#' Calculate Treatment Capacity Sensitivity Indicator
#'
#' Scores each water supply's treatment capacity from 0 (no or unknown
#' treatment) to 5 (full treatment train for its source type), using the
#' WSD `treatment_system_type` field. Adapted from Lynette Schwanger's
#' exploration (NPS_WSD/archive/exploration/NPS_WSD_Treatment_Indicator.Rmd).
#'
#' Steps:
#' 1. Each listed treatment is mapped to one or more treatment categories
#'    (disinfection, filtration, oxidation, sediment removal, specific
#'    contaminant removal, pH adjustment, corrosion control, emergency).
#' 2. Sources are scored on a groundwater or surface water track. GWUDI,
#'    spring, and rainwater sources are scored on both and keep the higher score.
#' 3. Coastal sources (within 50 km of the coast, from `coast_dist_m`) need
#'    reverse osmosis to score 4 or 5, because of saltwater intrusion risk.
#'
#' Groundwater track (surface water track in brackets):
#'   5 = oxidation [sediment removal] + filtration + disinfection + an add-on,
#'       and RO if coastal. Add-ons are sediment removal [oxidation], specific
#'       contaminant removal, pH adjustment, or corrosion control.
#'   4 = oxidation [sediment removal] + filtration + disinfection, and RO if coastal
#'   3 = filtration + disinfection (including a coastal full train without RO)
#'   2 = filtration or disinfection
#'   1 = only emergency, oxidation, sediment, contaminant, pH, or corrosion treatment
#'   0 = untreated or unknown treatment
#'
#' Unknown treatment (blank or "U" in the WSD) scores 0, an assumption of
#' minimum treatment capacity, and is flagged in `treatment_unknown`. Leaving
#' it NA would count as zero vulnerability in calc_vulnerability_index().
#'
#' Sources get an NA score, with the reason in `score_note`, when they are
#' inactive, non-potable, proposed, reclaimed, or treated by a municipal provider.
#' Only `active == "inactive"` (out of use for roughly 5-10+ years) counts as
#' inactive; "temporarily inactive" sources are scored like active ones.
#'
#' A higher raw score means more treatment capacity and lower sensitivity, so
#' norm_ and rank_ columns are inverted (higher = more sensitive).
#'
#' @param water_supplies sf or data frame of water supplies, usually
#'   `get_wsd(rapid = TRUE)`. Needs wsd_source_id, park_unit, park_name,
#'   classification, active, source_type, gwudi, and treatment_system_type.
#' @param coast_geometry file path to coastal_site_geometry.csv from
#'   process_coastal_geometry(). Sources with a non-NA `coast_dist_m` are coastal.
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with water supply ID, park info, treatment details, flags,
#'   and treatment capacity indicators (raw_, norm_, rank_)
#'
calc_treatment_capacity_sensitivity <- function(water_supplies,
                                                coast_geometry = here("data", "rapid", "inputs", "coast_comp", "coastal_site_geometry.csv"),
                                                save = FALSE,
                                                output_dir = NULL,
                                                verbose = TRUE) {

  required_cols <- c("wsd_source_id", "park_unit", "park_name", "classification", "active",
                     "source_type", "gwudi", "treatment_system_type")

  # Input validation
  stopifnot(
    "water_supplies is missing required WSD columns" = all(required_cols %in% names(water_supplies)),
    "coast_geometry file not found" = file.exists(coast_geometry),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )

  start_time <- Sys.time()

  treatment_category_map <- treatment_category_lookup()

  groundwater_types <- c("unknown groundwater", "deep groundwater", "shallow groundwater")
  surface_types     <- c("lake/reservoir", "river/stream", "ocean")
  both_track_types  <- c("spring", "rainwater")

  # Map a comma/semicolon-separated treatment string to unique categories
  map_treatments <- function(x) {
    if (is.na(x) || str_squish(x) %in% c("", "U")) return(NA_character_)
    treatments <- str_squish(unlist(strsplit(x, "[,;]")))
    treatments <- treatments[treatments != ""]
    categories <- unname(treatment_category_map[treatments])
    categories[is.na(categories)] <- "unmatched"
    paste(unique(str_squish(unlist(strsplit(categories, ",")))), collapse = ", ")
  }

  # Score one track. `primary` is the pretreatment step that completes the train:
  # oxidation for groundwater, sediment removal for surface water.
  score_track <- function(primary, addon, fi, di, ox, se, scr, ph, cc, em, ro, coastal) {
    ro_ok <- !coastal | ro
    case_when(
      primary & fi & di & addon & ro_ok    ~ 5,
      primary & fi & di & ro_ok            ~ 4,
      fi & di                              ~ 3,
      fi | di                              ~ 2,
      em | ox | se | scr | ph | cc         ~ 1,
      TRUE                                 ~ 0
    )
  }

  # 1. Prepare water supplies --------------------------------------------------
  coastal_sites <- readr::read_csv(coast_geometry, show_col_types = FALSE) |>
    select(wsd_source_id, coast_dist_m) |>
    distinct(wsd_source_id, .keep_all = TRUE)

  supplies <- water_supplies |>
    sf::st_drop_geometry() |>
    distinct(wsd_source_id, .keep_all = TRUE) |>
    select(all_of(required_cols)) |>
    left_join(coastal_sites, by = "wsd_source_id") |>
    mutate(coastal = !is.na(coast_dist_m))

  if (verbose) {
    message("Scoring treatment capacity for ", nrow(supplies), " water supplies...")
  }

  # 2. Map treatments to categories --------------------------------------------
  supplies <- supplies |>
    mutate(treatment_category = map_chr(treatment_system_type, map_treatments))

  unmatched <- supplies |>
    dplyr::filter(str_detect(coalesce(treatment_category, ""), "unmatched")) |>
    pull(treatment_system_type) |>
    strsplit("[,;]") |>
    unlist() |>
    str_squish() |>
    setdiff(names(treatment_category_map))
  if (length(unmatched) > 0) {
    warning("Treatment types not in the category map (ignored in scoring): ",
            paste(unique(unmatched), collapse = ", "))
  }

  # 3. Score -------------------------------------------------------------------
  results <- supplies |>
    mutate(
      cat = coalesce(treatment_category, ""),
      has_ox  = str_detect(cat, "oxidation"),
      has_fi  = str_detect(cat, "filtration"),
      has_di  = str_detect(cat, "disinfection"),
      has_se  = str_detect(cat, "sediment_removal"),
      has_scr = str_detect(cat, "specific_contaminant_removal"),
      has_ph  = str_detect(cat, "ph_adjustment"),
      has_cc  = str_detect(cat, "corrosion_control"),
      has_em  = str_detect(cat, "emergency"),
      has_reverse_osmosis = str_detect(coalesce(treatment_system_type, ""), "Reverse Osmosis"),
      treatment_unknown = is.na(treatment_category),

      scoring_track = case_when(
        coalesce(gwudi, "") == "yes" | source_type %in% both_track_types | is.na(source_type) ~ "Both (higher score)",
        source_type %in% groundwater_types ~ "Groundwater",
        source_type %in% surface_types ~ "Surface water",
        TRUE ~ NA_character_
      ),

      score_groundwater = score_track(
        primary = has_ox, addon = has_se | has_scr | has_ph | has_cc,
        fi = has_fi, di = has_di, ox = has_ox, se = has_se, scr = has_scr,
        ph = has_ph, cc = has_cc, em = has_em, ro = has_reverse_osmosis, coastal = coastal),
      score_surface = score_track(
        primary = has_se, addon = has_ox | has_scr | has_ph | has_cc,
        fi = has_fi, di = has_di, ox = has_ox, se = has_se, scr = has_scr,
        ph = has_ph, cc = has_cc, em = has_em, ro = has_reverse_osmosis, coastal = coastal),

      # Reason a source isn't scored. The first matching condition wins.
      score_note = case_when(
        # Only long-term closures; "temporarily inactive" sources are scored
        coalesce(active, "active") == "inactive" ~ "Inactive",
        classification == "Proposed System" ~ "Proposed",
        classification == "Non-Potable (NonPOT)" ~ "Non-potable",
        source_type == "reclaimed" ~ "Reclaimed",
        str_detect(cat, "municipal_treatment") ~ "Municipal treatment",
        is.na(scoring_track) ~ "Unrecognized source type",
        TRUE ~ NA_character_
      ),

      raw_treatment_capacity = case_when(
        !is.na(score_note) ~ NA_real_,
        treatment_unknown ~ 0,
        cat == "untreated" ~ 0,
        scoring_track == "Groundwater" ~ score_groundwater,
        scoring_track == "Surface water" ~ score_surface,
        TRUE ~ pmax(score_groundwater, score_surface)
      ),
      treatment_unknown = treatment_unknown & is.na(score_note)
    ) |>
    select(wsd_source_id, park_unit, park_name, classification, source_type, gwudi,
           treatment_system_type, treatment_category, scoring_track, coastal, coast_dist_m,
           has_reverse_osmosis, treatment_unknown, score_note, raw_treatment_capacity)

  # Calculate standardized indicators using helper function.
  # Higher treatment capacity = lower sensitivity, so invert.
  final_results <- results |>
    add_standardized_indicators("raw_treatment_capacity", invert = TRUE)

  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    filename <- file.path(output_dir, paste0("treatment_capacity_sensitivity_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)

    if (verbose) {
      message("Results saved to: ", filename)
    }
  }

  # Print timing
  if (verbose) {
    message(sum(!is.na(final_results$raw_treatment_capacity)), " sources scored, ",
            sum(final_results$treatment_unknown), " with unknown treatment scored 0, ",
            sum(!is.na(final_results$score_note)), " not scored")
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nTreatment capacity sensitivity processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }

  return(final_results)
}


#' Treatment Category Lookup
#'
#' Maps each WSD `treatment_system_type` entry to one or more treatment
#' categories (comma-separated). Used by calc_treatment_capacity_sensitivity().
#' Mapping from Lynette Schwanger's treatment indicator exploration.
#'
#' @return named character vector: names are WSD treatment types, values are categories
#'
treatment_category_lookup <- function() {
  c(
    "Activated Carbon"                 = "specific_contaminant_removal",
    "Aeration"                         = "oxidation",
    "Air Stripper"                     = "oxidation",
    "Anti-scalent Injection"           = "corrosion_control",
    "Aqua Mag Injection"               = "corrosion_control, specific_contaminant_removal",
    "Arsenic Removal"                  = "specific_contaminant_removal",
    "Bicarbonate Feed"                 = "ph_adjustment, corrosion_control",
    "Calcite System"                   = "ph_adjustment, corrosion_control",
    "Caustic Soda Injection"           = "ph_adjustment",
    "Chlorine Tablet Feeder"           = "disinfection",
    "Coagulation and Flocculation"     = "sediment_removal",
    "Corrosion Control"                = "corrosion_control",
    "Disinfection"                     = "disinfection",
    "Ferric Chloride"                  = "sediment_removal",
    "Filtration"                       = "filtration",
    "Gas Chlorine Disinfection"        = "disinfection",
    "Hand Filter"                      = "emergency",
    "Hydrogen Peroxide"                = "oxidation",
    "Hydrogen Sulfide Treatment"       = "specific_contaminant_removal",
    "Iodine"                           = "emergency",
    "Iron Removal"                     = "specific_contaminant_removal",
    "Iron and Manganese Removal"       = "specific_contaminant_removal",
    "Iron and Manganese Sequestration" = "specific_contaminant_removal",
    "Liquid Chlorine Disinfection"     = "disinfection",
    "Municipal Treatment"              = "municipal_treatment",
    "Nitrate Removal"                  = "specific_contaminant_removal",
    "Oxidation"                        = "oxidation",
    "Oxygen Injection"                 = "oxidation",
    "Ozone Treatment"                  = "oxidation",
    "Phosphate Injection"              = "corrosion_control",
    "Potassium Permanganate"           = "oxidation",
    "Pre-chlorination"                 = "oxidation",
    "Radon Removal"                    = "specific_contaminant_removal",
    "Reverse Osmosis"                  = "filtration, specific_contaminant_removal, oxidation",
    "Sand and Sediment Removal"        = "sediment_removal",
    "Soda Ash Injection"               = "ph_adjustment, corrosion_control",
    "Soda ash treatment"               = "ph_adjustment, corrosion_control",
    "Sulfur Removal"                   = "specific_contaminant_removal",
    "UV Disinfection"                  = "disinfection",
    "Untreated"                        = "untreated",
    "Water Softening"                  = "corrosion_control",
    "Zeolite Filtration"               = "filtration",
    "pH Adjustment"                    = "ph_adjustment"
  )
}
