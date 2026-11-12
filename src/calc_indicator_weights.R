#' @title Assign water supplies to indicator-weight source categories
#'
#' @description
#' Assigns each water supply to one of the six source categories used to look
#' up indicator weights: `"Streams & Rivers"`, `"Lakes & Reservoirs"`,
#' `"Springs & Shallow Groundwater"`, `"Deep Groundwater"`, `"Rainwater"` and
#' `"Ocean"`. Each category has its own weight for every indicator (for example,
#' flood indicators weigh streams fully and other sources partly), so this
#' table links each supply to the right set of weights.
#'
#' The category is assigned by these rules, in order:
#' - A groundwater or spring supply that is under the direct influence of
#'   surface water (`gwudi == "yes"`) is Springs & Shallow Groundwater,
#'   whatever its depth, because it behaves like a shallow source. The flag is
#'   ignored for other source types.
#' - `"river/stream"`, `"lake/reservoir"`, `"rainwater"` and `"ocean"` map
#'   directly to their category, and `"spring"` to Springs & Shallow
#'   Groundwater.
#' - `"shallow groundwater"`, `"deep groundwater"` and `"unknown groundwater"`
#'   are split by depth to water: shallower than 30.5 m (100 ft) is Springs &
#'   Shallow Groundwater, 30.5 m or deeper is Deep Groundwater. 100 ft is the
#'   start of the deepest depth-to-water class in the DRASTIC groundwater
#'   vulnerability rating (Aller et al. 1987).
#' - Anything else (for example a groundwater supply with no depth that isn't
#'   GWUDI) is `NA`.
#'
#' Depth to water for groundwater and spring supplies is, in order of
#' preference:
#' 1. the WSD's recorded `well_depth` (feet), assumed to sit 20% above the well
#'    bottom (`well_depth * 0.8`), converted to meters and capped at 200 m;
#' 2. otherwise, the modeled water table depth at the supply
#'    (`wtd_point_m` from [process_groundwater_depth()], from the CONUS-wide
#'    random-forest water table depth model, Zenodo record 18504963,
#'    <https://doi.org/10.1111/gwat.13362>).
#'
#' Other source types are given a depth of 0, which doesn't affect their
#' category.
#'
#' @references
#' Aller, L., Lehr, J. H., Petty, R., & Bennett, T. (1987). DRASTIC: A
#' standardized system to evaluate groundwater pollution potential using
#' hydrogeologic setting. Journal of the Geological Society of India, 29(1),
#' 23-37. <https://doi.org/10.17491/jgsi/1987/290112>
#'
#' @param water_supplies sf object. The water supplies, with `wsd_source_id`,
#'   `park_unit`, `park_name`, `source_type`, `well_depth` (recorded well
#'   depth, feet) and `gwudi` columns, e.g. from [npswavacore::get_wsd()]. A
#'   `wsd_source_id` that appears more than once is kept once. Required.
#' @param groundwater_data String or data frame. The modeled water table depth
#'   per supply, with `wsd_source_id` and `wtd_point_m` columns.
#'   * A path: read as the `wtd_all_sites.csv` written by
#'     [process_groundwater_depth()].
#'   * A data frame: used as given.
#'
#'   Required.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `indicator_weights_<date>.csv`, named with today's date,
#'     to `output_dir`, creating the folder if needed. Stops with an error if
#'     `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for the CSV file. Ignored
#'   when `save = FALSE`. Default: `NULL`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints a start message, the path of the saved file and the run
#'     time.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id` and
#'   `weight_category` (one of the six categories above, or `NA`).
#' @export
#' @seealso [process_groundwater_depth()], which makes `groundwater_data`;
#'   [calc_groundwater_depth_sensitivity()] for the depth-to-water sensitivity
#'   indicator.
#' @examplesIf dir.exists(npswavacore::wava_data_path("rapid", "inputs", "groundwater_depth"))
#' ws <- npswavacore::get_wsd(park = "ZION")
#' weights <- calc_indicator_weights(
#'   ws,
#'   groundwater_data = npswavacore::wava_data_path("rapid", "inputs",
#'                                                  "groundwater_depth",
#'                                                  "wtd_all_sites.csv")
#' )
#' dplyr::count(weights, weight_category)
calc_indicator_weights <- function(water_supplies,
                                   groundwater_data,
                                   save = FALSE,
                                   output_dir = NULL,
                                   verbose = TRUE) {
  
  # Input validation
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include source_type" = "source_type" %in% names(water_supplies),
    "water_supplies must include well_depth" = "well_depth" %in% names(water_supplies),
    "save is TRUE but output_dir not provided" = !save || !is.null(output_dir)
  )
  
  # Start timing
  if (verbose) {
    message("Starting categorical source type identification...")
    start_time <- Sys.time()
  }
  
  # Load pre-processed modeled WTD if a file path was given
  if (is.character(groundwater_data)) {
    groundwater_data <- read_csv(groundwater_data, show_col_types = FALSE)
  }
  
  # Score all water supplies. Groundwater/spring source types get a depth
  # from the well-depth/modeled-WTD priority chain below; all other
  # (non-groundwater/spring) source types are hardcoded to a depth of 0.
  
  all_supplies <- water_supplies %>%
    sf::st_drop_geometry() %>%
    dplyr::select(wsd_source_id, park_unit, park_name, source_type,
                  well_depth, gwudi) %>%
    dplyr::mutate(well_depth = as.numeric(well_depth)) %>%
    dplyr::distinct(wsd_source_id, .keep_all = TRUE)
  
  # Join modeled water table depth, then prefer the WSD's recorded well depth
  # (ground truth), assumed to sit 20% shallower than the well bottom, over
  # the modeled WTD.
  
  results <- all_supplies %>%
    dplyr::left_join(groundwater_data %>% dplyr::select(wsd_source_id, wtd_point_m),
                     by = "wsd_source_id") %>%
    dplyr::mutate(is_gw_spring = grepl("groundwater|spring", source_type),
                  well_depth_m = if_else((well_depth * 0.8 * 0.3048) >= 200, 200, (well_depth * 0.8 * 0.3048)),
                  raw_groundwater_depth_m = dplyr::case_when(!is_gw_spring ~ 0,
                                                             !is.na(well_depth_m) ~ well_depth_m,
                                                             !is.na(wtd_point_m) ~ wtd_point_m,
                                                             .default = NA_real_),
                  raw_groundwater_depth_source = dplyr::case_when(!is_gw_spring ~ "surface_zero",
                                                                  !is.na(well_depth_m) ~ "well_depth",
                                                                  !is.na(wtd_point_m) ~ "modeled_wtd",
                                                                  .default = NA_character_)) %>%
    
    # Bin the raw depth into a DRASTIC-style sensitivity score, then standardize
    # both the raw depth and the bin score separately
    
    mutate(groundwater_drastic = case_when(between(raw_groundwater_depth_m, 0, 1.5) ~ 1,
                                           between(raw_groundwater_depth_m, 1.5, 4.5) ~ .9,
                                           between(raw_groundwater_depth_m, 4.5, 9.1) ~ .7,
                                           between(raw_groundwater_depth_m, 9.1, 15.2) ~ .5,
                                           between(raw_groundwater_depth_m, 15.2, 22.9) ~ .3,
                                           between(raw_groundwater_depth_m, 22.9, 30.5) ~ .2,
                                           raw_groundwater_depth_m >= 30.5 ~ .1,
                                           .default = NA)) %>%
    # GWUDI groundwater and springs are checked first: groundwater under the
    # direct influence of surface water behaves like a shallow source whatever
    # its depth. is_gw_spring limits this to groundwater and spring source
    # types, so a surface-water supply flagged GWUDI keeps its own category.
    mutate(weight_category = case_when(is_gw_spring & gwudi == "yes" ~ "Springs & Shallow Groundwater",
                                       source_type == 'spring' ~ "Springs & Shallow Groundwater",
                                       source_type == 'river/stream' ~ "Streams & Rivers",
                                       source_type == 'rainwater' ~ "Rainwater",
                                       source_type == 'lake/reservoir' ~ "Lakes & Reservoirs",
                                       source_type == 'ocean' ~ "Ocean",
                                       source_type %in% c("shallow groundwater", "deep groundwater", 'unknown groundwater') & raw_groundwater_depth_m < 30.5 ~ "Springs & Shallow Groundwater",
                                       source_type %in% c("shallow groundwater", "deep groundwater", 'unknown groundwater') & raw_groundwater_depth_m >= 30.5 ~ "Deep Groundwater")) %>%
    select(wsd_source_id, weight_category)
  
  
  # Save results
  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    filename <- file.path(output_dir, paste0("indicator_weights_", Sys.Date(), ".csv"))
    write_csv(results, filename)
    
    if (verbose) {
      message("Results saved to: ", filename)
    }
  }
  
  # Print timing
  if (verbose) {
    end_time <- Sys.time()
    total_time <- end_time - start_time
    message("\nCategorical source type processing completed in ",
            round(total_time, 2), " ", attr(total_time, "units"))
  }
  
  return(results)
}
