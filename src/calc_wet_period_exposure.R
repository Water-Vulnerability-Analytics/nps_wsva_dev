#' @title Wet spell exposure: projected change in the longest wet spell
#'
#' @description
#' Calculates the change in the yearly maximum number of consecutive wet
#' days (CWD) at each water supply between a historical and a future period,
#' from the LOCA2 precipitation-extremes data processed by
#' [process_extreme_precip_data()]. For each model and ensemble member, the
#' historical mean comes from its historical run and the future mean from
#' each scenario run. The change is the future mean minus the historical
#' mean, in days and as a percentage of the historical mean.
#'
#' The ensemble is then summarized per site: the runs closest to the 10th,
#' 50th and 90th percentiles of each change, and the percentage of runs with
#' longer wet spells. Longer wet spells mean more exposure: more wet days in a
#' row saturate soils and raise flood and water-quality stress. Only the
#' percent change is standardized, so the days are not counted twice.
#'
#' @param water_supplies sf object. The water supplies, with `wsd_source_id`
#'   and `park_unit` columns, e.g. from [npswavacore::get_wsd()]. Only their IDs
#'   are used, to pick rows from the processed files. Required.
#' @param period_data String. The folder of parquet files from
#'   [process_extreme_precip_data()]. Stops with an error if it doesn't exist
#'   or the files lack the expected columns. Default:
#'   `rapid/inputs/loca_cmip6_precip_periods_water_supplies` in the shared
#'   data folder, see [npswavacore::wava_data_path()].
#' @param historical Character or numeric vector of 2 years, e.g.
#'   `c(1980, 2010)`. The first and last year of the historical period, read
#'   from the historical runs (1950-2014). Required.
#' @param future Character or numeric vector of 2 years, e.g.
#'   `c(2035, 2065)`. The first and last year of the future period, read from
#'   the scenario runs. Required.
#' @param save Logical. Whether to write the results to a CSV file.
#'   * `TRUE`: writes `wet_period_exposure_<date>.csv`, named with today's
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
#'   `park_unit` and:
#'   * `raw_cwd_change_p10`, `_p50`, `_p90`: the change in mean CWD (days) of
#'     the model run closest to that percentile;
#'   * `raw_cwd_percent_change_p10`, `_p50`, `_p90`: the same for the percent
#'     change;
#'   * `raw_cwd_change_models_showing_increase`: the percentage of runs with
#'     longer wet spells; `raw_cwd_change_model_count`: the number of runs;
#'   * `norm_*` and `rank_*` versions of the percent changes and of
#'     `raw_cwd_change_models_showing_increase`, from
#'     [add_standardized_indicators()]. Longer spells get a higher score.
#' @export
#' @seealso [calc_dry_period_exposure()], the paired dry spell indicator;
#'   [calc_wet_period_sensitivity()] for the historical trend;
#'   [process_extreme_precip_data()], which makes the input files.
#' @examplesIf dir.exists(npswavacore::wava_data_path("rapid", "inputs", "loca_cmip6_precip_periods_water_supplies"))
#' ws <- npswavacore::get_wsd(park = "ZION")
#' wet_exp <- calc_wet_period_exposure(ws,
#'                                     historical = c(1980, 2010),
#'                                     future = c(2035, 2065))
#' wet_exp[, c("wsd_source_id", "raw_cwd_percent_change_p50",
#'             "raw_cwd_change_models_showing_increase")]

calc_wet_period_exposure <- function(
    water_supplies,
    period_data = npswavacore::wava_data_path("rapid", "inputs", "loca_cmip6_precip_periods_water_supplies/"),
    historical,
    future,
    save = FALSE,
    output_dir = NULL,
    verbose = TRUE) {

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

  if (!dir.exists(period_data)) {
    stop("Processed precipitation-period directory not found: ", period_data)
  }

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
    message("Starting consecutive wet-day exposure calculation...")
    start_time <- Sys.time()
  }

  period_ds <- arrow::open_dataset(period_data, format = "parquet")

  required_cols <- c(
    "wsd_source_id",
    "park_unit",
    "year",
    "model",
    "ssp",
    "run",
    "wet_spell_max_length"
  )

  missing_cols <- setdiff(required_cols, names(period_ds))

  if (length(missing_cols)) {
    stop(
      "period_data is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  ws_ids <- unique(water_supplies$wsd_source_id)

  annual_data <- period_ds |>
    dplyr::select(dplyr::all_of(required_cols)) |>
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
    ) |>
    dplyr::collect() |>
    data.table::as.data.table()

  if (!nrow(annual_data)) {
    stop("No wet-period records remained after filtering.")
  }

  annual_data[
    !is.finite(wet_spell_max_length),
    wet_spell_max_length := NA_real_
  ]

  historical_means <- annual_data[
    ssp == "historical",
    .(
      cwd_historical = mean(wet_spell_max_length, na.rm = TRUE),
      hist_years_used = sum(!is.na(wet_spell_max_length))
    ),
    by = .(wsd_source_id, park_unit, model, run)
  ]

  future_means <- annual_data[
    ssp != "historical",
    .(
      cwd_future = mean(wet_spell_max_length, na.rm = TRUE),
      future_years_used = sum(!is.na(wet_spell_max_length))
    ),
    by = .(wsd_source_id, park_unit, model, ssp, run)
  ]

  wet_change <- merge(
    future_means,
    historical_means,
    by = c("wsd_source_id", "park_unit", "model", "run"),
    all = FALSE
  )

  wet_change[
    ,
    `:=`(
      cwd_change = cwd_future - cwd_historical,
      cwd_percent_change = data.table::fifelse(
        is.finite(cwd_historical) & cwd_historical > 0,
        100 * (cwd_future - cwd_historical) / cwd_historical,
        NA_real_
      )
    )
  ]

  wet_change <- wet_change[
    hist_years_used > 0L &
      future_years_used > 0L &
      is.finite(cwd_change) &
      is.finite(cwd_percent_change)
  ]

  if (!nrow(wet_change)) {
    stop(
      "No paired historical/future wet-period percent changes could be calculated."
    )
  }

  summary_dt <- wet_change[
    ,
    {
      change_quantiles <- stats::quantile(
        cwd_change,
        probs = c(0.10, 0.50, 0.90),
        na.rm = TRUE
      )
      percent_quantiles <- stats::quantile(
        cwd_percent_change,
        probs = c(0.10, 0.50, 0.90),
        na.rm = TRUE
      )

      closest_value <- function(values, target) {
        values[which.min(abs(values - target))]
      }

      .(
        raw_cwd_change_p10 =
          closest_value(cwd_change, change_quantiles[[1]]),
        raw_cwd_change_p50 =
          closest_value(cwd_change, change_quantiles[[2]]),
        raw_cwd_change_p90 =
          closest_value(cwd_change, change_quantiles[[3]]),
        raw_cwd_percent_change_p10 =
          closest_value(cwd_percent_change, percent_quantiles[[1]]),
        raw_cwd_percent_change_p50 =
          closest_value(cwd_percent_change, percent_quantiles[[2]]),
        raw_cwd_percent_change_p90 =
          closest_value(cwd_percent_change, percent_quantiles[[3]]),
        raw_cwd_change_models_showing_increase =
          100 * sum(cwd_percent_change > 0) / .N,
        raw_cwd_change_model_count = .N
      )
    },
    by = .(wsd_source_id, park_unit)
  ]

  final_results <- summary_dt |>
    # Longer projected wet spells indicate greater exposure. Standardize only
    # percent change so the parallel absolute-day fields are not double-counted.
    add_standardized_indicators("raw_cwd_percent_change_p") |>
    add_standardized_indicators(
      "raw_cwd_change_models_showing_increase"
    ) |>
    dplyr::select(wsd_source_id, park_unit, dplyr::everything())

  if (save) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    output_path <- file.path(
      output_dir,
      paste0("wet_period_exposure_", Sys.Date(), ".csv")
    )

    readr::write_csv(final_results, output_path)

    if (verbose) message("Results saved to: ", output_path)
  }

  if (verbose) {
    elapsed <- Sys.time() - start_time
    message(
      "Consecutive wet-day exposure completed in ",
      round(elapsed, 2), " ", attr(elapsed, "units")
    )
  }

  dplyr::as_tibble(final_results)
}
