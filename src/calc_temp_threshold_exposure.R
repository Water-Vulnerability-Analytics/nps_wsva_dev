#' @title Warm-water exposure: projected change in temperature thresholds
#'
#' @description
#' Calculates the change in three yearly temperature-threshold metrics at each
#' water supply between a historical and a future period, from the CMIP6-LOCA2
#' thresholds data processed by [process_temp_threshold_metrics()]. They are
#' air-temperature surrogates for surface water warming toward conditions that
#' favour harmful algal blooms, for the Water Quality temperature indicator:
#' - `su`: summer days (maximum above 25 C), hot-day frequency. Its threshold
#'   is near the growth optima reported for many bloom-forming cyanobacteria
#'   (Reinl et al. 2023; Joehnk et al. 2008).
#' - `gdd10`: growing degree days above 10 C, cumulative heat. Sustained
#'   warmth tracks water-column warming better than one hot day; degree days
#'   above 20 C best explained bloom timing across 580 Quebec lakes
#'   (Ratte-Fortin et al. 2020), so a 10 C base is low.
#' - `gsl`: growing season length (5 C threshold), warm-season length. It
#'   reaches about 365 days in the warmest places, where it can't show
#'   further warming.
#'
#' For each model, scenario and run, the historical and future means come from
#' the same file, which holds the historical run (1950-2014) and the scenario
#' (2015-2100) in one series. The change is the future mean minus the
#' historical mean, in the metric's units and as a percentage of the
#' historical mean.
#'
#' The ensemble is then summarized per site: the runs closest to the 10th,
#' 50th and 90th percentiles of the absolute change, and the percentage of
#' runs showing an increase. More warm days, heat or season mean more
#' exposure.
#'
#' Only the absolute change is standardized. Percent change inflates cool,
#' high-elevation supplies whose historical mean is near zero: a Rocky
#' Mountain supply with under 2 summer days a year shows over +1,000% for 26
#' more days, and for GDD10 percent change ranks supplies nearly opposite to
#' absolute change. Bloom thresholds are fixed biological limits, not relative
#' to the local climate, so the days or degree days added are what matter.
#' Percent change is reported for the same runs, and is `NA` where the
#' historical mean is zero.
#'
#' Air temperature stands in for surface water temperature, which lags and
#' damps it. These metrics apply to surface water supplies, and they don't
#' represent nutrient loading, the main control on bloom biomass, so they are
#' modifiers of bloom risk.
#'
#' @references
#' Joehnk, K. D., Huisman, J., Sharples, J., Sommeijer, B., Visser, P. M., &
#' Stroom, J. M. (2008). Summer heatwaves promote blooms of harmful
#' cyanobacteria. Global Change Biology, 14, 495-512.
#' <https://doi.org/10.1111/j.1365-2486.2007.01510.x>
#'
#' Ratte-Fortin, C., Chokmani, K., & Laurion, I. (2020). Spatiotemporal
#' variability in phytoplankton bloom phenology in eastern Canadian lakes
#' related to physiographic, morphologic, and climatic drivers. Environments,
#' 7, 77. <https://doi.org/10.3390/environments7100077>
#'
#' Reinl, K. L., et al. (2023). Blooms also like it cold. Limnology and
#' Oceanography Letters, 8, 546-564. <https://doi.org/10.1002/lol2.10316>
#'
#' @param water_supplies sf object. The water supplies, with `wsd_source_id`
#'   and `park_unit` columns, e.g. from [npswavacore::get_wsd()]. Only their IDs
#'   are used, to pick rows from the processed files. Required.
#' @param temp_threshold_data String. The parent folder of the `SU`, `GDD10`
#'   and `GSL` parquet folders from [process_temp_threshold_metrics()]. Stops
#'   with an error if it doesn't exist or a folder's files lack the expected
#'   columns; a missing metric folder gives a warning and that metric's columns
#'   are left out. Default:
#'   `rapid/inputs/loca_cmip6_temp_threshold_water_supplies` in the shared
#'   data folder, see [npswavacore::wava_data_path()].
#' @param historical Character or numeric vector of 2 years, e.g.
#'   `c(1980, 2010)`. The first and last year of the historical period. Gives
#'   a warning if it ends after 2014, since those years are scenario-forced.
#'   Required.
#' @param future Character or numeric vector of 2 years, e.g.
#'   `c(2035, 2065)`. The first and last year of the future period. Stops with
#'   an error if it starts before the historical period ends. Required.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `temp_threshold_exposure_<date>.csv`, named with today's
#'     date, to `output_dir`, creating the folder if needed. Stops with an
#'     error if `output_dir` is `NULL`.
#'   * `FALSE`: writes nothing.
#'
#'   Default: `FALSE`.
#' @param output_dir String or `NULL`. The folder for the CSV file.
#'   * A folder path: where the file is written when `save = TRUE`. Ignored
#'     when `save = FALSE`.
#'   * `NULL`: only allowed with `save = FALSE`.
#'
#'   Default: `NULL`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints each step, the run time and the path of the saved file.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return A tibble, one row per water supply, with `wsd_source_id`,
#'   `park_unit`, `historical_period`, `future_period` and, for each metric
#'   (`su`, `gdd10`, `gsl`):
#'   * `raw_<metric>_change_p10`, `_p50`, `_p90`: the change in the mean of
#'     the model run closest to that percentile;
#'   * `raw_<metric>_percent_change_p10`, `_p50`, `_p90`: the percent change
#'     of the same runs;
#'   * `raw_future_<metric>_p90`: the future mean of the p90 run;
#'   * `<metric>_historical`: the historical mean over all runs;
#'   * `raw_<metric>_models_showing_increase`: the percentage of runs with an
#'     increase; `<metric>_model_count`: the number of runs;
#'     `<metric>_members_with_percent_change`: the number with a percent
#'     change;
#'   * `norm_*` and `rank_*` versions of the absolute changes and of
#'     `raw_<metric>_models_showing_increase`, from
#'     [add_standardized_indicators()]. Larger increases get a higher score.
#' @export
#' @seealso [calc_temp_exposure()] for the change in mean temperature;
#'   [process_temp_threshold_metrics()], which makes the input files.
#' @examplesIf dir.exists(npswavacore::wava_data_path("rapid", "inputs", "loca_cmip6_temp_threshold_water_supplies"))
#' ws <- npswavacore::get_wsd(park = "ZION")
#' temp_threshold_exp <- calc_temp_threshold_exposure(ws,
#'                                                    historical = c(1980, 2010),
#'                                                    future = c(2035, 2065))
#' temp_threshold_exp[, c("wsd_source_id", "raw_su_change_p50",
#'                        "raw_gdd10_change_p50", "raw_gsl_change_p50")]

calc_temp_threshold_exposure <- function(
    water_supplies,
    temp_threshold_data = npswavacore::wava_data_path("rapid", "inputs", "loca_cmip6_temp_threshold_water_supplies/"),
    historical,
    future,
    save = FALSE,
    output_dir = NULL,
    verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # How this works, in brief. For each water supply and each climate model
  # run, we compare three warm-weather measures between a past period and a
  # future period: days above 25 °C (SU), growing degree days above 10 °C
  # (GDD10) and growing season length (GSL). Then we summarize across all the
  # model runs: a low (10th percentile), middle (median) and high (90th
  # percentile) projection, plus the share of runs that show an increase.
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

  if (!dir.exists(temp_threshold_data)) {
    stop("Processed temperature-threshold directory not found: ",
      temp_threshold_data)
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

  # Each file holds one model run from 1950 to 2100: past years first, then
  # the future scenario. Both periods are read from that same file, so if they
  # overlapped, the shared years would count on both sides of the change.
  if (future_start_year <= hist_end_year) {
    stop("The future period must start after the historical period ends.")
  }

  # From 2015 on, each file follows its own future scenario. A past period
  # that runs past 2014 would therefore differ from one scenario to the next.
  if (hist_end_year > 2014L) {
    warning(
      "The historical window ends after 2014, so it includes scenario-forced ",
      "years and differs by scenario."
    )
  }

  if (verbose) {
    message("Starting temperature-threshold exposure calculation...")
    start_time <- Sys.time()
  }

  # Step 3. The three measures (each also the name of its folder of processed
  # files), the short name each gets in the output columns, and which kind of
  # change is scored. All three are scored on the change in their own units
  # (days, or degree-days), not on percent change: some cold mountain supplies
  # have almost no days above 25 °C today, so even a small rise would look
  # like an enormous percentage.
  metric_spec <- list(
    list(variable = "SU", label = "su", standardize = "absolute"),
    list(variable = "GDD10", label = "gdd10", standardize = "absolute"),
    list(variable = "GSL", label = "gsl", standardize = "absolute")
  )

  # The columns that together identify one model run at one supply. The
  # scenario is part of it because each file, and so each past average,
  # belongs to one scenario.
  ws_ids <- unique(water_supplies$wsd_source_id)
  group_keys <- c("wsd_source_id", "park_unit", "model", "ssp", "run")

  # Step 4. Work through the three measures one at a time. Each has its own
  # folder of files, and a supply that is missing one measure still gets
  # results for the others.
  results <- NULL

  for (spec in metric_spec) {
    variable <- spec$variable
    label <- spec$label
    standardize <- spec$standardize
    metric_dir <- file.path(temp_threshold_data, variable)

    if (!dir.exists(metric_dir)) {
      warning(
        "No processed data for ", variable, " at ", metric_dir,
        "; its columns will be absent from the result."
      )
      next
    }

    if (verbose) message("  Summarizing ", variable)

    # Point at this measure's processed files without loading them yet, and
    # check they have the columns we need.
    metric_ds <- arrow::open_dataset(metric_dir, format = "parquet")

    required_cols <- c(
      "wsd_source_id", "park_unit", "year", "model", "ssp", "run", variable
    )
    missing_cols <- setdiff(required_cols, names(metric_ds))

    if (length(missing_cols)) {
      stop(
        metric_dir, " is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        ". Was it written by process_temp_threshold_metrics()?"
      )
    }

    # Load only our supplies and the years inside the two periods, and treat
    # impossible values (such as infinity) as missing.
    annual_data <- metric_ds |>
      dplyr::select(dplyr::all_of(required_cols)) |>
      dplyr::filter(
        wsd_source_id %in% ws_ids,
        (year >= hist_start_year & year <= hist_end_year) |
          (year >= future_start_year & year <= future_end_year)
      ) |>
      dplyr::collect() |>
      data.table::as.data.table()

    if (!nrow(annual_data)) {
      warning("No ", variable, " records remained after filtering.")
      next
    }

    annual_data[!is.finite(get(variable)), (variable) := NA_real_]

    # Average each model run over the past period and over the future period.
    historical_means <- annual_data[
      year >= hist_start_year & year <= hist_end_year,
      .(
        metric_historical = mean(get(variable), na.rm = TRUE),
        hist_years_used = sum(!is.na(get(variable)))
      ),
      by = group_keys
    ]

    future_means <- annual_data[
      year >= future_start_year & year <= future_end_year,
      .(
        metric_future = mean(get(variable), na.rm = TRUE),
        future_years_used = sum(!is.na(get(variable)))
      ),
      by = group_keys
    ]

    # Put the two averages side by side and work out the change: future minus
    # past, both in the measure's own units and as a percentage of the past.
    # The percentage is left blank when the past average is zero.
    metric_change <- merge(
      future_means,
      historical_means,
      by = group_keys,
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
    # fine: a supply with no days above 25 °C in the past is a real result,
    # not an error.
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
        # scored, and the percentage is reported for the same run so the two
        # can be read side by side.
        rank_values <- if (identical(standardize, "percent")) {
          metric_percent_change
        } else {
          metric_abs_change
        }
        # If ranking on percentages and no run has one (every past average was
        # zero), rank on the change in units instead of returning nothing.
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
          sum(is.finite(metric_percent_change))
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
    stop("No temperature-threshold metric produced a usable result.")
  }

  # Step 5. Rescale the scored columns so different indicators can be compared
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

    final_results <- final_results |>
      add_standardized_indicators(prefix) |>
      add_standardized_indicators(
        paste0("raw_", label, "_models_showing_increase")
      )
  }

  # Step 6. Record which periods were compared, then save to CSV if asked.
  final_results <- final_results |>
    dplyr::mutate(
      historical_period = paste0(hist_start_year, "-", hist_end_year),
      future_period = paste0(future_start_year, "-", future_end_year)
    ) |>
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())

  if (save) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    output_path <- file.path(
      output_dir,
      paste0("temp_threshold_exposure_", Sys.Date(), ".csv")
    )

    readr::write_csv(final_results, output_path)

    if (verbose) message("Results saved to: ", output_path)
  }

  if (verbose) {
    elapsed <- Sys.time() - start_time
    message(
      "Temperature-threshold exposure completed in ",
      round(elapsed, 2), " ", attr(elapsed, "units")
    )
  }

  dplyr::as_tibble(final_results)
}
