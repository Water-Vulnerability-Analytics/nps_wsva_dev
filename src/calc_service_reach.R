#' Calculate Water System Service Reach
#'
#' Assigns each water supply a service reach category (Wide, Moderate, Focused,
#' or Not assessed) with a reason. Reach describes how many people a system
#' serves and its role at the park; it doesn't rank systems as more or less
#' important. This isn't a scored indicator, so it isn't
#' normalized; it's used to filter or add context to the rapid results.
#'
#' Each system is rated on two things, then combined with park visitation:
#' - Scale: population served from EPA SDWA (EPA size breaks at 500 and 3,300).
#'   Municipal supplies get no population, because SDWA reports the whole utility.
#'   The WSD `source_type` decides whether the park runs the system:
#'   - Physical sources (groundwater, spring, surface water, rainwater, ...): park-run.
#'   - `consecutive`: the park buys piped water and runs its own storage and
#'     distribution under its own PWS. Park-run, with that PWS's population.
#'   - `hauled`: delivered into park storage, so park-run, unless `pws_number` is
#'     the supplier's PWS (its PWSID appears in `supplied_from`, or SDWA lists a
#'     local or state government owner). Water hauled from another park system
#'     (`supplied_from` gives a WSD system ID) is always park water.
#'   - `municipal`, `tribal`, `purchased/concessioner`: the park is a customer on
#'     the supplier's system. "Municipal/purchased", with no population.
#' - Park role: how the system compares with the park's other visitor-facing,
#'   park-run systems (Primary = largest, Secondary = 30%+ share, Minor, Staff only).
#'
#' Park roles and municipal counts need every active source at the park, including
#' municipal supplies that aren't on the rapid list, so they're calculated from
#' `all_water_supplies` before the output is limited to `water_supplies`.
#'
#' Activity status comes from the WSD only:
#' - Sources that are inactive in the WSD (out of use for roughly 5-10 years or
#'   more) are left out of every calculation. They're still returned, as "Not
#'   assessed" with the reason "Inactive in WSD", so every source in
#'   `water_supplies` has a row.
#' - Sources that are temporarily inactive in the WSD (recently or possibly
#'   temporarily closed) are kept and assessed like active sources.
#' - Sources that are active in the WSD but whose PWS isn't active in SDWA are
#'   kept and assessed like any other active source, using SDWA's last reported
#'   population. They're flagged in `review_sdwa_inactive` so the WSD and SDWA
#'   records can be reconciled.
#'
#' See Explore/explore_service_reach.Rmd for the logic, results, and the
#' items that need data clarification.
#'
#' @param water_supplies sf or data frame of the sources to return, usually
#'   `get_wsd(rapid = TRUE)`. Only `wsd_source_id` is used.
#' @param all_water_supplies all CONUS sources, including municipal, purchased,
#'   and hauled. Defaults to `get_wsd(conus = TRUE, fill_empty_coords = TRUE)`.
#' @param sdwa_file path to `sdwa_public_water_systems.csv`, the EPA SDWA public
#'   water system table saved by [process_sdwa_data()]. Defaults to
#'   `rapid/inputs` in the shared data folder. Its SDWA release
#'   (`SUBMISSIONYEARQUARTER`) and download date are added to the output as
#'   `sdwa_release` and `sdwa_download_date`.
#' @param visitation data frame of monthly visitation with UnitCode, Year, and
#'   RecreationVisitors, as returned by [npswavacore::get_unit_visitation()].
#'   Defaults to downloading `visitation_years` for every park with an active
#'   source, so it needs internet.
#' @param visitation_years years averaged for recreation visits. The default,
#'   2021-2025, leaves out 2020, when COVID closures lowered visitation.
#' @param tier_low_max,tier_moderate_max population breaks between Low, Moderate,
#'   and High scale.
#' @param assumed_pop_non_pws population given to non-public systems with no SDWA match.
#' @param assumed_pop_pws population given to PWS with no usable SDWA population.
#' @param secondary_share_min share of the park's visitor-facing population needed
#'   to be a Secondary system.
#' @param park_vis_high,park_vis_moderate annual recreation visits for a busy and a
#'   moderately busy park.
#' @param return_details logical. If TRUE, returns a list with the results, all
#'   CONUS sources (active and inactive), and the settings used.
#' @param save logical, whether to save output
#' @param output_dir character, directory path for saving files
#' @param verbose logical, whether to show progress messages
#'
#' @return tibble with one row per source in `water_supplies`,
#'   with `service_reach`, `service_reach_reason`, `wide_reach`, the
#'   supporting scale, role, and visitation columns, and `review_*` flags.
#'   A list if `return_details = TRUE`.
#'
calc_service_reach <- function(water_supplies,
                               all_water_supplies = NULL,
                               sdwa_file = npswavacore::wava_data_path("rapid", "inputs",
                                                                       "sdwa_public_water_systems.csv"),
                               visitation = NULL,
                               visitation_years = 2021:2025,
                               tier_low_max = 500,
                               tier_moderate_max = 3300,
                               assumed_pop_non_pws = 24,
                               assumed_pop_pws = 25,
                               secondary_share_min = 0.3,
                               park_vis_high = 1000000,
                               park_vis_moderate = 100000,
                               return_details = FALSE,
                               save = FALSE,
                               output_dir = NULL,
                               verbose = TRUE) {

  # Input validation
  stopifnot(
    "water_supplies must contain wsd_source_id" = "wsd_source_id" %in% names(water_supplies),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )

  start_time <- Sys.time()

  pws_categories <- c("Low", "Moderate", "High")
  park_run_categories <- c("Non-PWS", pws_categories)
  # Source types where the park is a customer on the supplier's system
  customer_pattern <- regex("municipal|tribal|purchased", ignore_case = TRUE)
  # Source types where the water comes from outside the park
  bought_pattern <- regex("municipal|tribal|purchased|consecutive|hauled", ignore_case = TRUE)
  # supplied_from values that point to another WSD system (e.g. "HALE_WS_001")
  wsd_system_pattern <- regex("\\b[A-Z]{4}_WS_\\d{3}\\b")

  # Keyword groups for the free-text point_of_use field. The first match wins, so a
  # source serving both visitors and staff counts as visitor-facing. Only
  # "Staff/admin" changes the result: those systems are left out of the role comparison.
  pou_patterns <- c(
    "Lodging"         = "lodg|hotel|inn\\b|cabin|motel|resort|village",
    "Campground"      = "camp|\\brv\\b",
    "Visitor day use" = "visitor|\\bvc\\b|picnic|restroom|comfort|beach|bath|restaurant|store|marina|museum|day use|trailhead|fountain|concession|shower|snack|guest|rafter",
    "Parkwide"        = "all park|park facilities|entire park|headquarters|\\bhq\\b",
    "Staff/admin"     = "resid|hous|employee|staff|quarters|dorm|maint|office|ranger station|shop"
  )

  classify_point_of_use <- function(x) {
    type <- rep(NA_character_, length(x))
    for (name in names(pou_patterns)) {
      match <- str_detect(coalesce(x, ""), regex(pou_patterns[[name]], ignore_case = TRUE))
      type <- if_else(is.na(type) & match, name, type)
    }
    case_when(!is.na(type) ~ type,
              is.na(x) | trimws(x) == "" ~ "Missing",
              TRUE ~ "Other")
  }

  # WSD codes that differ from the visitation codes (same as calc_visitation_sensitivity)
  visitation_code_crosswalk <- c("ANAC" = "NCPE", "NACE" = "NCPE", "NAMA" = "NCPC", "VALR" = "PERL")

  # 1. Import all sources and set aside those inactive in the WSD ---------------
  if (verbose) message("Reading water supplies and SDWA data...")

  if (is.null(all_water_supplies)) {
    all_water_supplies <- get_wsd(conus = TRUE, fill_empty_coords = TRUE)
  }
  if (!file.exists(sdwa_file)) {
    stop("SDWA file not found: ", sdwa_file, "\nRun process_sdwa_data() first ",
         "(see rapid_data_processing.Rmd).", call. = FALSE)
  }
  sdwa <- readr::read_csv(sdwa_file,
                          col_types = readr::cols(.default = readr::col_character(),
                                                  POPULATION_SERVED_COUNT = readr::col_double()),
                          show_col_types = FALSE)
  # Which SDWA release the results use; EPA replaces its "latest" file each quarter
  sdwa_release <- paste(sort(unique(sdwa$SUBMISSIONYEARQUARTER)), collapse = ", ")
  sdwa_download_date <- paste(sort(unique(sdwa$download_date)), collapse = ", ")

  all_sources <- all_water_supplies |>
    sf::st_drop_geometry() |>
    distinct(wsd_source_id, .keep_all = TRUE) |>
    select(park_unit, park_name, water_system_name, wsd_system_id, wsd_source_id, classification,
           system_type, source_type, use_category, active, pws_number, supplied_from, point_of_use)

  # Activity status comes from the WSD. Only sources marked "inactive" (out of use
  # for roughly 5-10+ years) are left out of every calculation below and added back
  # at the end as "Not assessed". "temporarily inactive" sources are recent or
  # possibly temporary closures, so they stay in. A missing status counts as active.
  wsd_inactive <- all_sources |>
    dplyr::filter(coalesce(active, "active") == "inactive")

  sources <- all_sources |>
    dplyr::filter(coalesce(active, "active") != "inactive")

  sdwa_clean <- sdwa |>
    mutate(
      pws_number_clean = toupper(trimws(PWSID)),
      pws_type_desc = case_when(
        PWS_TYPE_CODE == "CWS"    ~ "Community Water System",
        PWS_TYPE_CODE == "NTNCWS" ~ "Non-Transient Non-Community",
        PWS_TYPE_CODE == "TNCWS"  ~ "Transient Non-Community",
        TRUE ~ "Unknown/Other"
      ),
      owner_type_desc = case_when(
        OWNER_TYPE_CODE == "F" ~ "Federal Government",
        OWNER_TYPE_CODE == "L" ~ "Local Government",
        OWNER_TYPE_CODE == "S" ~ "State Government",
        OWNER_TYPE_CODE == "P" ~ "Private",
        OWNER_TYPE_CODE == "M" ~ "Mixed Ownership",
        OWNER_TYPE_CODE == "N" ~ "Native American",
        TRUE ~ "Other"
      )
    ) |>
    select(pws_number_clean,
           pws_name            = PWS_NAME,
           pws_activity_status = PWS_ACTIVITY_CODE,
           pws_type_desc,
           owner_type_desc,
           population_served   = POPULATION_SERVED_COUNT)

  # One SDWA record per source. A few sources list two PWS numbers, so prefer the
  # record that's active in SDWA, then the largest population served.
  sdwa_by_source <- sources |>
    select(wsd_source_id, pws_number) |>
    separate_rows(pws_number, sep = ",\\s*") |>
    mutate(pws_number_clean = toupper(trimws(pws_number))) |>
    inner_join(sdwa_clean, by = "pws_number_clean") |>
    arrange(wsd_source_id, desc(pws_activity_status == "A"), desc(population_served)) |>
    distinct(wsd_source_id, .keep_all = TRUE) |>
    select(-pws_number)

  # SDWA activity status doesn't remove anything. A source that's active in the WSD
  # stays in, even if its PWS is inactive in SDWA; it uses SDWA's last reported
  # population and is flagged in review_sdwa_inactive (step 6).

  # Mean annual recreation visits. Total visits also count non-recreation traffic
  # such as commuters, which inflates parkways.
  if (is.null(visitation)) {
    visitation_codes <- unique(coalesce(visitation_code_crosswalk[sources$park_unit], sources$park_unit))
    if (verbose) message("Downloading visitation for ", length(visitation_codes), " park units...")
    # The API takes a comma-separated list of units, so request them in batches
    visitation <- split(visitation_codes, ceiling(seq_along(visitation_codes) / 50)) |>
      purrr::map(\(units) npswavacore::get_unit_visitation(units = units,
                                                         startYear = min(visitation_years),
                                                         endYear = max(visitation_years))) |>
      purrr::list_rbind()
  }

  park_visits <- visitation |>
    dplyr::filter(Year %in% visitation_years) |>
    group_by(UnitCode, Year) |>
    summarize(annual_visits = sum(RecreationVisitors, na.rm = TRUE), .groups = "drop") |>
    group_by(UnitCode) |>
    summarize(park_rec_visits = mean(annual_visits))

  # 2-3. Scale category and point of use ---------------------------------------
  if (verbose) message("Assigning scale and park role to ", nrow(sources), " active sources...")

  # The first matching condition wins.
  system_categories <- sources |>
    left_join(sdwa_by_source, by = "wsd_source_id") |>
    mutate(
      has_sdwa_pop = !is.na(population_served) & population_served > 0,
      # Water from another park system (supplied_from gives its WSD ID) is park water
      in_park_transfer = str_detect(coalesce(supplied_from, ""), wsd_system_pattern),
      # Water bought or hauled from outside the park (consecutive, hauled, municipal,
      # tribal, or purchased source)
      outside_source = str_detect(coalesce(source_type, ""), bought_pattern) & !in_park_transfer,
      # Hauled water goes into park storage, but pws_number can be the supplier's PWS.
      # It is when its PWSID appears in supplied_from or a town or state owns it.
      hauled_supplier_pws = str_detect(coalesce(source_type, ""), regex("hauled", ignore_case = TRUE)) &
        !in_park_transfer & !is.na(pws_number_clean) &
        (str_detect(toupper(coalesce(supplied_from, "")), fixed(coalesce(pws_number_clean, "@"))) |
           coalesce(owner_type_desc, "") %in% c("Local Government", "State Government")),
      # The park is a customer on the supplier's system: SDWA reports the whole
      # utility, not the park's portion, so these get no population. Consecutive
      # systems are park-run and keep their own PWS's population.
      is_municipal = (str_detect(coalesce(source_type, ""), customer_pattern) & !in_park_transfer) |
        hauled_supplier_pws,
      population_estimate = case_when(
        classification %in% c("Non-Potable (NonPOT)", "Proposed System") | is_municipal ~ NA_real_,
        has_sdwa_pop ~ population_served,
        classification == "Non-Public Water System (NP)" ~ assumed_pop_non_pws,
        TRUE ~ assumed_pop_pws # most likely a missing or incorrect PWS number
      ),
      population_source = case_when(
        is.na(population_estimate) ~ NA_character_,
        has_sdwa_pop ~ "SDWA",
        TRUE ~ "Assumed"
      ),
      system_category = case_when(
        classification == "Proposed System" ~ "Proposed",
        classification == "Non-Potable (NonPOT)" ~ "Non-potable",
        is_municipal ~ "Municipal/purchased",
        !has_sdwa_pop & classification == "Non-Public Water System (NP)" ~ "Non-PWS",
        population_estimate <= tier_low_max ~ "Low",
        population_estimate <= tier_moderate_max ~ "Moderate",
        TRUE ~ "High"
      ),
      pou_type = classify_point_of_use(point_of_use),
      # Each PWS counts once at a park, however many sources it has. Systems without
      # a matched population fall back to the WSD system ID.
      unit_id = if_else(has_sdwa_pop, pws_number_clean, wsd_system_id),
      visitation_code = coalesce(visitation_code_crosswalk[park_unit], park_unit),
      # Classification defines non-potable; flag systems where use_category disagrees
      review_potable_conflict = !is.na(use_category) & (
        (system_category == "Non-potable" & str_detect(use_category, "human-drinking")) |
          (system_category %in% park_run_categories & !str_detect(use_category, "human-drinking"))
      ),
      # Active in the WSD, but the matched PWS isn't active in SDWA. Kept and assessed;
      # flagged so the WSD and SDWA records can be reconciled.
      review_sdwa_inactive = !is.na(pws_activity_status) & pws_activity_status != "A",
      # SDWA population may not describe the park: park-run systems buying or hauling
      # water can serve others, and non-federal CWS listed as park-run may be utilities
      review_population = system_category %in% pws_categories & (
        outside_source |
          (coalesce(pws_type_desc, "") == "Community Water System" &
             coalesce(owner_type_desc, "") != "Federal Government")
      )
    )

  # 4. Park role ----------------------------------------------------------------
  # Each park's own potable systems are compared with each other. Staff-only systems
  # are left out, so the role reflects the systems visitors rely on.
  park_roles <- system_categories |>
    dplyr::filter(system_category %in% park_run_categories) |>
    group_by(park_unit, unit_id) |>
    summarize(population_estimate = first(population_estimate),
              visitor_facing = any(pou_type != "Staff/admin"),
              .groups = "drop") |>
    group_by(park_unit) |>
    mutate(n_park_run = n()) |>
    group_by(park_unit, visitor_facing) |>
    mutate(share_visitor_facing = if_else(visitor_facing, population_estimate / sum(population_estimate), NA_real_),
           is_largest = visitor_facing & population_estimate == max(population_estimate),
           n_tied_primary = sum(is_largest)) |>
    ungroup() |>
    mutate(park_role = case_when(
      !visitor_facing ~ "Staff only",
      is_largest ~ "Primary",
      share_visitor_facing >= secondary_share_min ~ "Secondary",
      TRUE ~ "Minor"
    )) |>
    select(park_unit, unit_id, n_park_run, share_visitor_facing, n_tied_primary, park_role)

  # Active municipal connections per park. This counts connections, not volume supplied.
  park_municipal <- system_categories |>
    dplyr::filter(system_category == "Municipal/purchased") |>
    distinct(park_unit, wsd_system_id) |>
    count(park_unit, name = "n_municipal")

  # 5-6. Park context and service reach ------------------------------------------
  all_results <- system_categories |>
    left_join(park_roles, by = c("park_unit", "unit_id")) |>
    left_join(park_municipal, by = "park_unit") |>
    left_join(park_visits, by = c("visitation_code" = "UnitCode")) |>
    group_by(park_unit) |>
    mutate(n_park_run = if (all(is.na(n_park_run))) 0L else max(n_park_run, na.rm = TRUE)) |>
    ungroup() |>
    mutate(
      n_municipal = replace_na(n_municipal, 0L),
      park_role = factor(coalesce(park_role, "Not park-run"),
                         levels = c("Primary", "Secondary", "Minor", "Staff only", "Not park-run")),
      park_visitation_class = case_when(
        is.na(park_rec_visits) ~ "No data",
        park_rec_visits >= park_vis_high ~ "High",
        park_rec_visits >= park_vis_moderate ~ "Moderate",
        TRUE ~ "Low"
      ),
      # Only PWS can qualify through the park role: non-PWS systems serve fewer than 25
      # people per day by definition. Parks with as many municipal connections as
      # park-run systems are left out, because their own systems may not be their main supply.
      role_counts = system_category %in% pws_categories & n_municipal < n_park_run,
      # The first matching condition wins.
      service_reach_reason = case_when(
        system_category == "Non-potable" ~ "Non-potable",
        system_category == "Municipal/purchased" ~ "Municipal or purchased",
        system_category == "Proposed" ~ "Proposed",
        system_category == "High" ~ "High scale",
        role_counts & park_role == "Primary" & park_visitation_class == "High" ~ "Primary system at a busy park",
        system_category == "Moderate" ~ "Moderate scale",
        role_counts & park_role == "Primary" & park_visitation_class == "Moderate" ~ "Primary system at a moderately busy park",
        role_counts & park_role == "Secondary" & park_visitation_class == "High" ~ "Secondary system at a busy park",
        park_role == "Staff only" ~ "Staff-only system",
        system_category == "Non-PWS" ~ "Non-public system",
        TRUE ~ "Small system"
      ),
      service_reach = case_when(
        service_reach_reason %in% c("Non-potable", "Municipal or purchased", "Proposed") ~ "Not assessed",
        service_reach_reason %in% c("High scale", "Primary system at a busy park") ~ "Wide",
        service_reach_reason %in% c("Moderate scale", "Primary system at a moderately busy park",
                                 "Secondary system at a busy park") ~ "Moderate",
        TRUE ~ "Focused"
      ) |>
        factor(levels = c("Wide", "Moderate", "Focused", "Not assessed")),
      wide_reach = service_reach == "Wide",
      # Review flags
      review_staff_only_pws = park_role == "Staff only" & system_category %in% pws_categories,
      review_tied_primary = park_role == "Primary" & coalesce(n_tied_primary, 1L) > 1,
      review_assumed_pws = population_source == "Assumed" & system_category %in% pws_categories,
      review_population_vs_visits = population_source == "SDWA" & !is.na(park_rec_visits) &
        population_estimate >= park_rec_visits / 365,
      across(starts_with("review_"), \(x) coalesce(x, FALSE)),
      system_category = factor(system_category, levels = c("High", "Moderate", "Low", "Non-PWS",
                                                           "Municipal/purchased", "Non-potable", "Proposed"))
    ) |>
    select(wsd_source_id, park_unit, park_name, water_system_name, wsd_system_id,
           service_reach, service_reach_reason, wide_reach,
           system_category, population_estimate, population_source, park_role,
           share_visitor_facing, n_tied_primary, n_park_run, n_municipal,
           park_rec_visits, park_visitation_class,
           classification, source_type, use_category, point_of_use, pou_type, unit_id,
           outside_source, in_park_transfer, supplied_from,
           pws_number_clean, pws_name, pws_type_desc, owner_type_desc, pws_activity_status,
           population_served, starts_with("review_"))

  # 7. Add back WSD-inactive sources and limit to the requested sources ----------
  # Inactive sources weren't assessed, but they're returned so every source in
  # water_supplies has a service reach.
  all_results <- all_results |>
    bind_rows(
      wsd_inactive |>
        transmute(wsd_source_id, park_unit, park_name, water_system_name, wsd_system_id,
                  service_reach = factor("Not assessed", levels = levels(all_results$service_reach)),
                  service_reach_reason = "Inactive in WSD",
                  wide_reach = FALSE,
                  classification, source_type, use_category, point_of_use)
    ) |>
    left_join(select(all_sources, wsd_source_id, active), by = "wsd_source_id") |>
    relocate(active, .after = service_reach_reason) |>
    mutate(across(starts_with("review_"), \(x) coalesce(x, FALSE)))

  # Roles and municipal counts above used all active sources, so they still
  # reflect each park's full set of systems.
  requested_ids <- unique(water_supplies$wsd_source_id)
  final_results <- all_results |>
    dplyr::filter(wsd_source_id %in% requested_ids) |>
    mutate(sdwa_release = sdwa_release, sdwa_download_date = sdwa_download_date)

  missing_ids <- setdiff(requested_ids, final_results$wsd_source_id)
  if (length(missing_ids) > 0) {
    warning(length(missing_ids), " sources in water_supplies aren't in all_water_supplies ",
            "and have no service reach: ", paste(head(missing_ids, 10), collapse = ", "))
  }

  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    filename <- file.path(output_dir, paste0("service_reach_", Sys.Date(), ".csv"))
    write_csv(final_results, filename)

    if (verbose) {
      message("Results saved to: ", filename)
    }
  }

  # Print timing
  if (verbose) {
    message(nrow(final_results), " sources returned (",
            sum(final_results$service_reach_reason == "Inactive in WSD"), " inactive in the WSD, ",
            sum(final_results$review_sdwa_inactive), " active in the WSD but not in SDWA); ",
            sum(final_results$wide_reach), " Wide; SDWA release ", sdwa_release,
            " (downloaded ", sdwa_download_date, ")")
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nService reach processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }

  if (return_details) {
    return(list(
      results = final_results,
      all_sources = all_results,
      settings = list(
        sdwa_release = sdwa_release, sdwa_download_date = sdwa_download_date,
        tier_low_max = tier_low_max, tier_moderate_max = tier_moderate_max,
        assumed_pop_non_pws = assumed_pop_non_pws, assumed_pop_pws = assumed_pop_pws,
        secondary_share_min = secondary_share_min, park_vis_high = park_vis_high,
        park_vis_moderate = park_vis_moderate, visitation_years = visitation_years,
        pou_patterns = pou_patterns
      )
    ))
  }

  return(final_results)
}
