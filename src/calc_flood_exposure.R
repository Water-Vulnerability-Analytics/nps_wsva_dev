#' Calculate Flood Exposure Indicators from LOCA2 Precipitation Extremes
#'
#' Calculates projected change in three flood-relevant precipitation extremes
#' from the processed revised LOCA2 (v20240915) precipitation-extremes dataset.
#' For each model/SSP/run, change is the difference between historical mean and
#' future mean; the model value closest to the p10, p50 and p90 of the local
#' ensemble is retained, following the same approach as
#' `calc_dry_period_exposure()` and the other water-balance exposure
#' indicators.
#'
#' Three metrics are returned because flood-producing rain differs with
#' catchment size and storm structure, and no single daily index captures all
#' of it. The set follows Crespi et al. (2020, section 4.2.2), whose European
#' hazard-index framework uses daily-precipitation indices for the
#' extreme-precipitation hazard. Two studies that test
#' ETCCDI indices against observed floods support the two depth metrics
#' (Avila et al. 2016; Pesce et al. 2022).
#'
#' - **`prmax1day`** (single-day intensity). Annual highest 1-day total, ETCCDI
#'   Rx1day, in mm. 
#' - **`prmax5day`** (multi-day accumulation). Annual highest 5-day total,
#'   ETCCDI Rx5day, in mm.
#' - **`pr_days_above_nonzero_99th`** (heavy-rain frequency). Annual days over
#'   the model's own 1951-1980 wet-day 99th percentile, in days. Matches the
#'   frequency-of-extreme-precipitation index of Crespi et al. (2020, index
#'   15; medium priority), which they list as applicable to water management.
#'   It measures repeated heavy-rain stress independent of absolute volume.
#'
#' Larger positive changes indicate greater exposure, so no indicator is
#' inverted.
#'
#' **Scope.** Precipitation is used here as a direct indicator of pluvial and
#' fluvial *meteorological* exposure. This is a macro-scale screening comparison
#' of relative vulnerability shifts across water-supply points, not an
#' asset-level engineering or insurance assessment. It does not model hydraulic
#' routing, local topography, channel geometry or built drainage capacity, so it
#' identifies supplies under increasing rainfall stress but cannot determine
#' site-level flood depth or ponding. 
#' 
#' @references
#' Avila, A., Justino, F., Wilson, A., Bromwich, D., & Amorim, M. (2016).
#' Recent precipitation trends, flash floods and landslides in southern Brazil.
#' Environmental Research Letters, 11, 114029.
#' <https://doi.org/10.1088/1748-9326/11/11/114029>
#'
#' Crespi, A., Terzi, S., Cocuccioni, S., Zebisch, M., Berckmans, J., &
#' Fuessel, H.-M. (2020). Climate-related hazard indices for Europe. ETC/CCA
#' Technical Paper 1/2020, European Topic Centre on Climate Change Impacts,
#' Vulnerability and Adaptation.
#' <https://doi.org/10.25424/cmcc/climate_related_hazard_indices_europe_2020>
#'
#' Pesce, M., von Hardenberg, J., Claps, P., & Viglione, A. (2022).
#' Correlation between climate and flood indices in Northwestern Italy at
#' different temporal scales. Journal of Hydrology and Hydromechanics, 70(2),
#' 178-194. <https://doi.org/10.2478/johh-2022-0009>
#'
#' @param water_supplies An `sf` object containing `wsd_source_id` and
#'   `park_unit`.
#' @param flood_data Character. Directory containing parquet files created by
#'   `process_extreme_precip_data()`.
#' @param historical Character or numeric vector of length 2 giving the
#'   historical start/end years.
#' @param future Character or numeric vector of length 2 giving the future
#'   start/end years.
#' @param save Logical. If `TRUE`, writes results to disk.
#' @param output_dir Character. Output directory when `save = TRUE`.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#'
#' @return A tibble with one row per water supply and, for each of the three
#'   metrics, the p10/p50/p90 absolute and percent change, the share of models
#'   showing an increase, the model count, and the historical mean.
calc_flood_exposure <- function(
    water_supplies,
    flood_data = npswavacore::wava_data_path("rapid", "inputs", "loca_cmip6_precip_periods_water_supplies/"),
    historical,
    future,
    save = FALSE,
    output_dir = NULL,
    verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # How this works, in brief. For each water supply and each climate model
  # run, we compare three heavy-rain measures between a past period and a
  # future period. Then we summarize across all the model runs: a low (10th
  # percentile), middle (median) and high (90th percentile) projection, plus
  # the share of runs that show an increase.
  # ---------------------------------------------------------------------------
  #
  # Step 1. Check the inputs.
  stopifnot(
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" =
      "wsd_source_id" %in% names(water_supplies),
    "water_supplies must include park_unit" =
      "park_unit" %in% names(water_supplies),
    "historical must be length 2" = length(historical) == 2,
    "future must be length 2" = length(future) == 2,
    "save is TRUE but output_dir not provided" =
      !save || !is.null(output_dir)
  )

  if (!dir.exists(flood_data)) {
    stop("Processed flood-metrics directory not found: ", flood_data)
  }

  # Step 2. Turn the periods into start and end years. A year ("1980") or a
  # date ("1980-01-01") both work, because only the first four characters are
  # used.
  hist_start_year <- as.integer(substr(historical[1], 1, 4))
  hist_end_year <- as.integer(substr(historical[2], 1, 4))
  future_start_year <- as.integer(substr(future[1], 1, 4))
  future_end_year <- as.integer(substr(future[2], 1, 4))

  if (anyNA(c(
    hist_start_year,
    hist_end_year,
    future_start_year,
    future_end_year
  ))) {
    stop("Could not parse years from historical/future inputs.")
  }

  if (
    hist_start_year > hist_end_year ||
      future_start_year > future_end_year
  ) {
    stop("Period start years must be less than or equal to end years.")
  }

  if (verbose) {
    message("Starting flood exposure calculation...")
    start_time <- Sys.time()
  }

  # Step 3. The three measures, the short name each gets in the output
  # columns, and which kind of change is scored. The two rainfall amounts
  # (wettest day and wettest 5 days) are scored on percent change, so wet and
  # dry places can be compared fairly. The count of very heavy rain days is
  # scored on the change in days, because its past average is often close to
  # zero, and a percent change from near zero is meaningless.
  metric_spec <- list(
    list(variable = "prmax1day",
         label = "prmax1day",
         standardize = "percent"),
    list(variable = "pr_days_above_nonzero_99th",
         label = "pr99days",
         standardize = "absolute"),
    list(variable = "prmax5day",
         label = "prmax5day",
         standardize = "percent")
  )

  # Step 4. Point at the processed files without loading them yet, and check
  # they have every column we need.
  flood_ds <- arrow::open_dataset(flood_data, format = "parquet")

  metric_variables <- vapply(metric_spec, `[[`, character(1), "variable")

  required_cols <- c(
    "wsd_source_id",
    "park_unit",
    "year",
    "model",
    "ssp",
    "run",
    metric_variables
  )

  missing_cols <- setdiff(required_cols, names(flood_ds))

  if (length(missing_cols)) {
    stop(
      "flood_data is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      ". Was it written by process_extreme_precip_data()?"
    )
  }

  # Step 5. Load only the rows we need: our supplies, past years from each
  # model's historical run, and future years from its scenario runs.
  ws_ids <- unique(water_supplies$wsd_source_id)

  annual_data <- flood_ds %>%
    dplyr::select(dplyr::all_of(required_cols)) %>%
    dplyr::filter(
      wsd_source_id %in% ws_ids,
      (
        ssp == "historical" &
          year >= hist_start_year &
          year <= hist_end_year
      ) |
        (
          ssp != "historical" &
            year >= future_start_year &
            year <= future_end_year
        )
    ) %>%
    dplyr::collect() %>%
    data.table::as.data.table()

  if (!nrow(annual_data)) {
    stop("No flood-metric records remained after filtering.")
  }

  # Treat impossible values (such as infinity) as missing.
  for (variable in metric_variables) {
    annual_data[
      !is.finite(get(variable)),
      (variable) := NA_real_
    ]
  }

  # Step 6. Work through the three measures one at a time, so a supply that is
  # missing one measure still gets results for the other two.
  results <- NULL

  for (spec in metric_spec) {
    variable <- spec$variable
    label <- spec$label
    standardize <- spec$standardize

    if (verbose) message("  Summarizing ", variable)

    # Average each model run over the past period. Each model has one
    # historical run, which all of its future scenarios share.
    historical_means <- annual_data[
      ssp == "historical",
      .(
        metric_historical = mean(get(variable), na.rm = TRUE),
        hist_years_used = sum(!is.na(get(variable)))
      ),
      by = .(wsd_source_id, park_unit, model, run)
    ]

    # Average each model run over the future period, separately for each
    # scenario.
    future_means <- annual_data[
      ssp != "historical",
      .(
        metric_future = mean(get(variable), na.rm = TRUE),
        future_years_used = sum(!is.na(get(variable)))
      ),
      by = .(wsd_source_id, park_unit, model, ssp, run)
    ]

    # Pair each future run with the historical run of the same model, then
    # work out the change: future minus past, both in the measure's own units
    # and as a percentage of the past. The percentage is left blank when the
    # past average is zero.
    metric_change <- merge(
      future_means,
      historical_means,
      by = c("wsd_source_id", "park_unit", "model", "run"),
      all = FALSE
    )

    metric_change[
      ,
      `:=`(
        metric_abs_change = metric_future - metric_historical,
        metric_percent_change = data.table::fifelse(
          is.finite(metric_historical) & metric_historical > 0,
          100 * (metric_future - metric_historical) / metric_historical,
          NA_real_
        )
      )
    ]

    # Keep a run only if it had data in both periods. A blank percentage is
    # fine: for the heavy-rain day count, a past average of zero is a real
    # result, not an error.
    metric_change <- metric_change[
      hist_years_used > 0L &
        future_years_used > 0L &
        is.finite(metric_abs_change)
    ]

    if (!nrow(metric_change)) {
      warning(
        "No paired historical/future changes could be calculated for ",
        variable, "; its columns will be absent from the result."
      )
      next
    }

    summary_dt <- metric_change[
      ,
      {
        # Summarize the model runs for this supply. Rather than report an in-
        # between value that no model actually produced, we find the real run
        # closest to the low (10th percentile), middle (median) and high (90th
        # percentile) projection. Runs are ranked on the change that gets
        # scored (percent for the rainfall amounts, days for the day count),
        # and the other kind of change is reported for the same run so the two
        # can be read side by side.
        rank_values <- if (identical(standardize, "percent")) {
          metric_percent_change
        } else {
          metric_abs_change
        }
        # If no run has a percentage (every past average was zero), rank on
        # the change in units instead of returning nothing.
        if (!any(is.finite(rank_values))) rank_values <- metric_abs_change

        change_quantiles <- stats::quantile(
          rank_values,
          probs = c(0.10, 0.50, 0.90),
          na.rm = TRUE
        )

        # Find the run whose change is closest to a target value. Runs with a
        # blank value are skipped.
        closest_index <- function(values, target) {
          which.min(abs(values - target))
        }

        p10_i <- closest_index(rank_values, change_quantiles[[1]])
        p50_i <- closest_index(rank_values, change_quantiles[[2]])
        p90_i <- closest_index(rank_values, change_quantiles[[3]])

        n_percent <- sum(is.finite(metric_percent_change))

        # The values kept for each supply; their column names are set just
        # below.
        out <- list(
          metric_abs_change[p10_i],
          metric_abs_change[p50_i],
          metric_abs_change[p90_i],
          metric_percent_change[p10_i],
          metric_percent_change[p50_i],
          metric_percent_change[p90_i],
          metric_future[p90_i],
          mean(metric_historical, na.rm = TRUE),
          100 * sum(metric_abs_change > 0) / .N,
          .N,
          n_percent
        )

        names(out) <- c(
          paste0("raw_", label, "_change_p10"),
          paste0("raw_", label, "_change_p50"),
          paste0("raw_", label, "_change_p90"),
          paste0("raw_", label, "_percent_change_p10"),
          paste0("raw_", label, "_percent_change_p50"),
          paste0("raw_", label, "_percent_change_p90"),
          paste0("raw_future_", label, "_p90"),
          paste0(label, "_historical"),
          paste0("raw_", label, "_models_showing_increase"),
          paste0(label, "_model_count"),
          paste0(label, "_members_with_percent_change")
        )

        out
      },
      by = .(wsd_source_id, park_unit)
    ]

    # Add this measure's columns to the results so far.
    results <- if (is.null(results)) {
      summary_dt
    } else {
      merge(
        results,
        summary_dt,
        by = c("wsd_source_id", "park_unit"),
        all = TRUE
      )
    }
  }

  if (is.null(results)) {
    stop("No flood metric produced a usable result.")
  }

  # Step 7. Rescale the scored columns so different indicators can be compared
  # and combined. For each, add_standardized_indicators() adds a norm_ column
  # (0 for the lowest supply, 1 for the highest) and a rank_ column (the
  # supply's percentile, 0-100). Only the scored change and the share of runs
  # showing an increase are rescaled, so nothing is counted twice when
  # indicators are combined.
  final_results <- dplyr::as_tibble(results)

  for (spec in metric_spec) {
    label <- spec$label
    prefix <- if (identical(spec$standardize, "percent")) {
      paste0("raw_", label, "_percent_change_p")
    } else {
      paste0("raw_", label, "_change_p")
    }

    if (!any(startsWith(names(final_results), prefix))) next

    final_results <- final_results %>%
      add_standardized_indicators(prefix) %>%
      add_standardized_indicators(
        paste0("raw_", label, "_models_showing_increase")
      )
  }

  # Step 8. Record which periods were compared, then save to CSV if asked.
  final_results <- final_results %>%
    dplyr::mutate(
      historical_period = paste0(hist_start_year, "-", hist_end_year),
      future_period = paste0(future_start_year, "-", future_end_year)
    ) %>%
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())

  if (save) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    output_path <- file.path(
      output_dir,
      paste0("flood_exposure_", Sys.Date(), ".csv")
    )

    readr::write_csv(final_results, output_path)

    if (verbose) message("Results saved to: ", output_path)
  }

  if (verbose) {
    elapsed <- Sys.time() - start_time
    message(
      "Flood exposure completed in ",
      round(elapsed, 2), " ", attr(elapsed, "units")
    )
  }

  dplyr::as_tibble(final_results)
}
