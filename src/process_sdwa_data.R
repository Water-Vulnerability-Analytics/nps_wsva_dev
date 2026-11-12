#' @title EPA SDWA public water system attributes for the rapid workflow
#'
#' @description
#' Prepares the EPA Safe Drinking Water Act (SDWA) public water system table
#' that [calc_service_reach()] reads. Downloads EPA ECHO's
#' `SDWA_latest_downloads.zip` if it isn't already in `dest_dir` (or when
#' `overwrite = TRUE`), unzips it, and saves the columns the rapid workflow
#' uses from `SDWA_PUB_WATER_SYSTEMS.csv` to
#' `sdwa_public_water_systems.csv` in `output_dir`.
#'
#' EPA replaces the "latest" file each quarter, so the saved table keeps
#' `SUBMISSIONYEARQUARTER` and adds `download_date`, which record the SDWA
#' release a run used.
#'
#' A download is checked before it replaces anything: it goes to a temporary
#' file, the HTTP status must be a success, and the file must be a zip that
#' contains `SDWA_PUB_WATER_SYSTEMS.csv`. Only then does it replace the zip and
#' the unzipped folder. A failed download stops with an error and leaves the
#' existing files as they were.
#'
#' @param dest_dir String. The folder for the raw zip
#'   (`SDWA_latest_downloads.zip`) and its unzipped folder
#'   (`SDWA_latest_downloads/`). Created if it doesn't exist. Default:
#'   `all/Water_Supply` in the shared data folder, see
#'   [npswavacore::wava_data_path()].
#' @param output_dir String or `NULL`. The folder for
#'   `sdwa_public_water_systems.csv`. Default: `rapid/inputs` in the shared
#'   data folder, which is where [calc_service_reach()] looks for it. Only
#'   `NULL` with `save = FALSE`.
#' @param overwrite Logical. Whether to download a new copy of the zip.
#'   * `FALSE`: uses the zip already in `dest_dir`, and downloads only if
#'     there isn't one.
#'   * `TRUE`: downloads EPA's current release, replacing the zip and the
#'     unzipped folder once the download has been checked.
#'
#'   Default: `FALSE`.
#' @param url String. Where to download the zip. Default: EPA ECHO's SDWA
#'   download.
#' @param save Logical. Whether to write `sdwa_public_water_systems.csv` to
#'   `output_dir`, replacing any existing file. Default: `TRUE`.
#' @param verbose Logical. Whether to print progress. Default: `TRUE`.
#'
#' @return A tibble, one row per public water system, with
#'   `SUBMISSIONYEARQUARTER`, `PWSID`, `PWS_NAME`, `PWS_TYPE_CODE`,
#'   `OWNER_TYPE_CODE`, `PRIMARY_SOURCE_CODE`, `POPULATION_SERVED_COUNT`,
#'   `PWS_ACTIVITY_CODE`, `PWS_DEACTIVATION_DATE` and `download_date` (the
#'   date the zip was downloaded, from its file time).
#'
#' @seealso \url{https://echo.epa.gov/tools/data-downloads/sdwa-download-summary}
#'   for the SDWA dataset and field definitions.
process_sdwa_data <- function(
    dest_dir   = npswavacore::wava_data_path("all", "Water_Supply"),
    output_dir = npswavacore::wava_data_path("rapid", "inputs"),
    overwrite  = FALSE,
    url        = "https://echo.epa.gov/files/echodownloads/SDWA_latest_downloads.zip",
    save       = TRUE,
    verbose    = TRUE
) {

  stopifnot("save is TRUE but output_dir not provided" = !save || !is.null(output_dir))

  start_time <- Sys.time()
  pws_csv <- "SDWA_PUB_WATER_SYSTEMS.csv"
  keep_cols <- c("SUBMISSIONYEARQUARTER", "PWSID", "PWS_NAME", "PWS_TYPE_CODE",
                 "OWNER_TYPE_CODE", "PRIMARY_SOURCE_CODE", "POPULATION_SERVED_COUNT",
                 "PWS_ACTIVITY_CODE", "PWS_DEACTIVATION_DATE")

  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  zip_path  <- file.path(dest_dir, "SDWA_latest_downloads.zip")
  unzip_dir <- file.path(dest_dir, "SDWA_latest_downloads")
  pws_path  <- file.path(unzip_dir, pws_csv)

  # 1. Download, only when there's no zip or overwrite = TRUE ------------------
  downloaded <- FALSE
  if (!file.exists(zip_path) || overwrite) {
    if (verbose) message("Downloading the SDWA dataset (over 500 MB, this may take a while)...")
    tmp_zip <- tempfile(fileext = ".zip")
    on.exit(unlink(tmp_zip), add = TRUE)
    response <- httr::GET(url, httr::write_disk(tmp_zip, overwrite = TRUE),
                          if (verbose) httr::progress())
    # A failed request still writes the response body, so check before using it
    httr::stop_for_status(response, task = "download the SDWA dataset")
    contents <- tryCatch(utils::unzip(tmp_zip, list = TRUE)$Name, error = function(e) character(0))
    if (!pws_csv %in% basename(contents)) {
      stop("The download isn't a zip containing ", pws_csv, ". The existing files were left as they were.",
           call. = FALSE)
    }
    file.copy(tmp_zip, zip_path, overwrite = TRUE)
    downloaded <- TRUE
    if (verbose) message("Download complete: ", zip_path)
  } else if (verbose) {
    message("Using the existing zip (", format(file.mtime(zip_path), "%Y-%m-%d"),
            "). Use `overwrite = TRUE` to download EPA's current release.")
  }

  # 2. Unzip after a download, or when the folder is missing or incomplete -----
  if (downloaded || !file.exists(pws_path)) {
    if (verbose) message("Unzipping...")
    tmp_dir <- tempfile("sdwa_")
    on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
    utils::unzip(zip_path, exdir = tmp_dir)
    extracted <- list.files(tmp_dir, pattern = paste0("^", pws_csv, "$"), recursive = TRUE, full.names = TRUE)
    if (length(extracted) == 0) {
      stop(zip_path, " doesn't contain ", pws_csv, ". Run with `overwrite = TRUE` to download it again.",
           call. = FALSE)
    }
    # Replace the old folder only once the new files are extracted
    unlink(unzip_dir, recursive = TRUE)
    dir.create(unzip_dir)
    file.copy(list.files(dirname(extracted[1]), full.names = TRUE), unzip_dir, overwrite = TRUE)
    if (verbose) message("Unzipped to: ", unzip_dir)
  }

  # 3. Keep the columns the rapid workflow uses --------------------------------
  if (verbose) message("Reading: ", pws_path)
  sdwa <- readr::read_csv(pws_path, col_select = dplyr::all_of(keep_cols),
                          col_types = readr::cols(.default = readr::col_character(),
                                                  POPULATION_SERVED_COUNT = readr::col_double()),
                          show_col_types = FALSE) |>
    dplyr::mutate(download_date = as.Date(file.mtime(zip_path)))

  releases <- unique(sdwa$SUBMISSIONYEARQUARTER)
  if (verbose) message(nrow(sdwa), " public water systems, SDWA release ", paste(releases, collapse = ", "),
                       ", downloaded ", format(unique(sdwa$download_date)))

  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    out_path <- file.path(output_dir, "sdwa_public_water_systems.csv")
    readr::write_csv(sdwa, out_path)
    if (verbose) message("Results saved to: ", out_path)
  }

  if (verbose) message("SDWA processing completed in ",
                       round(difftime(Sys.time(), start_time, units = "secs"), 2), " secs")
  sdwa
}
