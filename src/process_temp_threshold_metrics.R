#' @title Extract LOCA2 temperature-threshold metrics at water supplies
#'
#' @description
#' Reads the CMIP6-LOCA2 thresholds collection (Alder 2024,
#' <https://doi.org/10.5066/P13OV6GY>) and writes the yearly metric at each
#' water supply to parquet files that [calc_temp_threshold_exposure()] reads.
#' `path` can be either:
#' - one `.tar.gz` archive as downloaded (one per metric, holding a NetCDF file
#'   for every model, scenario and run; 72 files and about 13 GB for SU, GDD10
#'   and GSL together). Files that already have a parquet are skipped without
#'   being unpacked; the rest are unpacked one batch at a time to a local
#'   temporary folder (not the synced shared data folder), processed, and
#'   deleted. This is the normal way to run it.
#' - one NetCDF file already unpacked from an archive.
#'
#' Each NetCDF file holds one metric for one model, scenario and run, with the
#' historical run (1950-2014) and the scenario (2015-2100) in one series. The
#' metric is one of:
#' - `SU`: summer days, the year's number of days with a maximum temperature
#'   above 25 C (ETCCDI SU), days;
#' - `GDD10`: growing degree days above 10 C, C d (stored as scaled integers;
#'   `ncdf4` applies the scale factor);
#' - `GSL`: growing season length (ETCCDI GSL, 5 C threshold), days.
#'
#' The model, scenario and run are read from each file's `model_name`,
#' `experiment` and `ripf` global attributes, since the file names
#' (`CMIP6-LOCA2_Thresholds_<METRIC>_<MODEL>.<SSP>.<RUN>_1950-2100_16thdeg_grid.nc`)
#' don't follow the precipitation collection's pattern. The extraction itself
#' is [process_extreme_precip_data()]: the mean over the 3 x 3 block of grid
#' cells around the supply's nearest cell, ignoring missing cells.
#'
#' Each parquet file is `<nc file name>.parquet` in `<output_dir>/<metric>`, one
#' folder per metric so each folder has one set of columns: `wsd_source_id`,
#' `park_unit`, `year`, `model`, `ssp`, `run`, `nc_file` and the metric.
#'
#' @param path String. Path to one
#'   `CMIP6-LOCA2_Thresholds_AllModels_grid_<METRIC>.tar.gz` archive, or to one
#'   NetCDF file unpacked from it, in `all/loca_cmip6_ext_temperature_products`
#'   in the shared data folder. Stops with an error if it doesn't exist, isn't
#'   a `.tar.gz` or `.nc` file, or holds a metric not in `metrics`. Required.
#' @param water_supplies sf object. The water supplies, as `POINT` geometries
#'   with `wsd_source_id` and `park_unit` columns, e.g. from
#'   [npswavacore::get_wsd()]. Stops with an error if any point is outside the
#'   LOCA2 grid. Required.
#' @param metrics Character vector. The metric names to accept. Default:
#'   `c("SU", "GDD10", "GSL")`.
#' @param grab_nearest,nearest_max_km Passed to
#'   [process_extreme_precip_data()]: whether to use the nearest cell with
#'   data, within `nearest_max_km` km, for supplies whose 3 x 3 block has no
#'   data. Default: `FALSE` and `200`.
#' @param output_dir String. The parent folder of the per-metric folders,
#'   created if needed. Default:
#'   `rapid/inputs/loca_cmip6_temp_threshold_water_supplies` in the shared
#'   data folder, see [npswavacore::wava_data_path()].
#' @param extract_dir String. A local folder for unpacked files when `path` is
#'   an archive, deleted on exit. It needs room for the archive's contents
#'   (about 9.3 GB for GDD10). Ignored for a NetCDF file. Default: a new
#'   folder in [tempdir()].
#' @param overwrite Logical. What to do with files already processed.
#'   * `TRUE`: processes every file again and replaces its parquet.
#'   * `FALSE`: skips files whose parquet exists (for an archive, without
#'     unpacking them).
#'
#'   Default: `FALSE`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints the archive, how many files are unpacked, each file and
#'     the run time.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return The parquet path(s), invisibly: one for a NetCDF file, one per
#'   NetCDF file in an archive.
#' @export
#' @seealso [calc_temp_threshold_exposure()], which reads the parquet files;
#'   [process_extreme_precip_data()], which does the extraction.
#' @examplesIf interactive() && dir.exists(npswavacore::wava_data_path("all", "loca_cmip6_ext_temperature_products"))
#' ws <- npswavacore::get_wsd(park = "ZION")
#' # A whole archive (the GSL archive is the smallest, 1.5 GB)
#' archive <- npswavacore::wava_data_path("all", "loca_cmip6_ext_temperature_products",
#'                                        "CMIP6-LOCA2_Thresholds_AllModels_grid_GSL.tar.gz")
#' process_temp_threshold_metrics(archive, ws, output_dir = tempdir())
#'
#' # Or one file unpacked from an archive
#' archive <- npswavacore::wava_data_path("all", "loca_cmip6_ext_temperature_products",
#'                                        "CMIP6-LOCA2_Thresholds_AllModels_grid_SU.tar.gz")
#' member <- "CMIP6-LOCA2_Thresholds_AllModels_grid/SU/CMIP6-LOCA2_Thresholds_SU_GFDL-CM4.ssp585.r1i1p1f1_1950-2100_16thdeg_grid.nc"
#' utils::untar(archive, files = member, exdir = tempdir())
#' out_file <- process_temp_threshold_metrics(file.path(tempdir(), member), ws,
#'                                            output_dir = tempdir())
#' head(arrow::read_parquet(out_file))

process_temp_threshold_metrics <- function(
    path,
    water_supplies,
    metrics = c("SU", "GDD10", "GSL"),
    grab_nearest = FALSE,
    nearest_max_km = 200,
    output_dir = npswavacore::wava_data_path("rapid", "inputs", "loca_cmip6_temp_threshold_water_supplies/"),
    extract_dir = tempfile("loca2_temp_thresholds_"),
    overwrite = FALSE,
    verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # How this works, in brief. The temperature measures (summer days, growing
  # degree days and growing season length) come as large compressed archives,
  # one per measure, each holding one NetCDF file per climate model run. Each
  # NetCDF file is laid out like the extreme-precipitation files, so the
  # actual extraction is handed to process_extreme_precip_data(). This
  # function handles the rest: unpacking archives a few files at a time,
  # working out which measure and which model run each file is, and saving
  # each measure in its own folder.
  # ---------------------------------------------------------------------------

  # Step 1. Check the inputs, and work out whether we were given an archive or
  # a single NetCDF file.
  stopifnot(
    "path must identify one existing file" = length(path) == 1 &&
      file.exists(path),
    "metrics must be a non-empty character vector" =
      is.character(metrics) && length(metrics) > 0
  )

  is_archive <- grepl("\\.(tar\\.gz|tgz)$", path, ignore.case = TRUE)
  is_netcdf <- grepl("\\.nc$", path, ignore.case = TRUE)

  if (!is_archive && !is_netcdf) {
    stop("path must be a .tar.gz archive or a .nc file: ", basename(path))
  }

  # Step 2. The work for one NetCDF file, used for a single file and for each
  # file unpacked from an archive.
  process_one_file <- function(nc_path) {
    # Look inside the file: which measure does it hold, and what do its
    # descriptive labels (global attributes) say?
    nc <- ncdf4::nc_open(nc_path)
    global_attributes <- ncdf4::ncatt_get(nc, 0)
    metric_variable <- intersect(names(nc$var), metrics)
    ncdf4::nc_close(nc)

    # Each file should hold exactly one of the measures we know about.
    if (length(metric_variable) != 1L) {
      stop(
        "Expected exactly one of ", paste(metrics, collapse = ", "),
        " in ", basename(nc_path), "; found ",
        if (length(metric_variable)) paste(metric_variable, collapse = ", ") else "none",
        "."
      )
    }

    # These file names don't reliably say which model, scenario and run they
    # are, so read that from the file's own labels, and stop if any label is
    # missing.
    required_attributes <- c("model_name", "experiment", "ripf")
    missing_attributes <- required_attributes[
      !vapply(
        required_attributes,
        function(x) {
          value <- global_attributes[[x]]
          is.character(value) && length(value) == 1L && nzchar(value)
        },
        logical(1)
      )
    ]

    if (length(missing_attributes)) {
      stop(
        basename(nc_path), " is missing global attributes: ",
        paste(missing_attributes, collapse = ", ")
      )
    }

    # Hand the file to the shared extractor. Each measure is saved in its own
    # folder (for example .../SU/), so every folder holds files with the same
    # columns.
    process_extreme_precip_data(
      nc_path = nc_path,
      water_supplies = water_supplies,
      variables = metric_variable,
      model_metadata = list(
        model = global_attributes$model_name,
        ssp = global_attributes$experiment,
        run = global_attributes$ripf
      ),
      grab_nearest = grab_nearest,
      nearest_max_km = nearest_max_km,
      output_dir = file.path(output_dir, metric_variable),
      overwrite = overwrite,
      verbose = verbose
    )
  }

  # Step 3. A single NetCDF file: process it and we're done.
  if (is_netcdf) {
    return(invisible(process_one_file(path)))
  }

  # Step 4. An archive. Unpacking a whole archive at once can need about 9 GB
  # of disk, so list what is inside, skip files already processed, unpack only
  # what is still needed, and delete each file as soon as it is done.
  if (verbose) message("Listing ", basename(path))

  # List the files in the archive without unpacking them.
  members <- utils::untar(path, list = TRUE)
  members <- members[grepl("\\.nc$", members)]

  if (!length(members)) {
    stop("No NetCDF files found in ", basename(path))
  }

  # Each file name says which measure it holds; check they are all ones we
  # expect.
  member_metric <- sub(
    "^CMIP6-LOCA2_Thresholds_([^_]+)_.*$",
    "\\1",
    basename(members)
  )
  unknown_metric <- !member_metric %in% metrics

  if (any(unknown_metric)) {
    stop(
      "Archive files do not match an accepted metric: ",
      paste(utils::head(basename(members[unknown_metric]), 3), collapse = ", ")
    )
  }

  # Where each file's result will be saved, so we can tell which are already
  # done.
  parquet_paths <- file.path(
    output_dir,
    member_metric,
    paste0(basename(members), ".parquet")
  )

  to_process <- if (overwrite) {
    rep(TRUE, length(members))
  } else {
    !file.exists(parquet_paths)
  }

  if (!any(to_process)) {
    if (verbose) {
      message("All ", length(members), " files already processed; skipping.")
    }
    return(invisible(parquet_paths))
  }

  # Unpack only the files still to do into a temporary folder, which is
  # deleted when the function finishes.
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(extract_dir, recursive = TRUE), add = TRUE)

  if (verbose) {
    message(
      "Extracting ", sum(to_process), " of ", length(members),
      " files to ", extract_dir
    )
  }

  extract_status <- utils::untar(
    path,
    files = members[to_process],
    exdir = extract_dir
  )

  if (!identical(as.integer(extract_status), 0L)) {
    stop("Extraction failed for ", basename(path))
  }

  # Process each unpacked file in turn, deleting it once it is done to free
  # disk space.
  for (member in members[to_process]) {
    process_one_file(file.path(extract_dir, member))
    unlink(file.path(extract_dir, member))
  }

  invisible(parquet_paths)
}
