#' @title Extract LOCA2 / MWBM values at water supplies from a NetCDF file
#'
#' @description
#' Reads one LOCA2 or MWBM (monthly water balance model) NetCDF file and
#' writes the values around each water supply to a parquet file. Run it once
#' per file; the `calc_*()` climate indicators read the folder of parquet
#' files.
#'
#' For each water supply it finds the nearest grid cell and keeps the 3 x 3
#' block of cells around it (fewer at the grid edge). It keeps the time steps
#' in the `historical` and `future` windows and converts units:
#' * precipitation stored as a flux (not `"mm"`) becomes monthly totals in mm
#'   (multiplied by 86400 and the days in the month);
#' * temperatures in kelvin become degrees C;
#' * daily Livneh precipitation is summed to monthly totals.
#'
#' The output file is `<nc file name>.parquet` in `output_dir`
#' (`<nc file name>_Tmin.parquet` or `_Tmax.parquet` for the "lusu"
#' temperature files), with one row per water supply, grid cell and time step
#' and columns `wsd_source_id`, `park_unit`, `time_units`, `nc_file`,
#' `cell_lat`, `cell_lon`, `date` and `value`. An existing file is replaced.
#'
#' @param nc_path String. Path to one NetCDF file, usually from
#'   `all/loca2_mwbm` in the shared data folder. Required.
#' @param water_supplies sf object. The water supplies, as points with
#'   `wsd_source_id` and `park_unit` columns, e.g. from [npswavacore::get_wsd()].
#'   They are reprojected to the file's CRS if needed. Required.
#' @param var String. The variable to extract: `"aet"`, `"deficit"`, `"pet"`,
#'   `"runoff"`, `"snow"`, `"swe"`, `"stor"`, `"pr"`, `"pr_tavg"`, `"tas"`,
#'   `"tasmin_tavg"` or `"tasmax_tavg"`. Some are read under another name
#'   inside the file:
#'   * `"pr"` is read as `pr_tavg`, except in the Jay Alder ("alder") files;
#'     `pr_tavg` is read as `PRCP` in the `livneh_unsplit_precip` files.
#'   * `"tasmin_tavg"` and `"tasmax_tavg"` are read as `Tmin` and `Tmax` in
#'     the "lusu" files.
#'   * `"swe"` is read as `snow`.
#'
#'   Stops with an error if `var` isn't supported or the variable isn't in
#'   the file. Required.
#' @param historical Character vector of 2 dates, `"YYYY-MM-DD"`. The start
#'   and end of the historical window. Only applied when the file's time
#'   units are "days since ..."; otherwise every time step is kept. Default:
#'   `c("1980-01-01", "2024-12-31")`.
#' @param future Character vector of 2 dates, `"YYYY-MM-DD"`. The start and
#'   end of the future window, applied the same way as `historical`. If no
#'   time step falls in either window, the function writes nothing and
#'   returns `NULL` with a message. Default: `c("2035-01-01", "2065-12-31")`.
#' @param output_dir String. The folder the parquet file is written to,
#'   created if it doesn't exist. Required.
#' @param grab_nearest Logical. What to do when none of a site's 3 x 3 cells
#'   has data (e.g. a coastal site over water).
#'   * `TRUE`: moves the block to the nearest cell with data within
#'     `nearest_max_km` and uses its 3 x 3 neighborhood instead. Sites with no
#'     data within that distance keep `NA`.
#'   * `FALSE`: keeps the site with `NA` values.
#'
#'   Default: `FALSE`.
#' @param nearest_max_km Numeric. The farthest, in km, `grab_nearest` looks for
#'   a cell with data. Ignored when `grab_nearest = FALSE`. Default: `200`,
#'   which reaches Dry Tortugas, about 180 km from the nearest LOCA2 cell with
#'   data.
#'
#' @return `NULL`, invisibly. Called for the parquet file it writes.
#' @export
#' @seealso [calc_runoff_exposure()], [calc_precip_exposure()],
#'   [calc_swe_exposure()], [calc_temp_exposure()], [process_drought_data()]
#'   and [calc_loca2_change_summary()], which read the parquet files.
#' @examplesIf interactive() && dir.exists(npswavacore::wava_data_path("all", "loca2_mwbm"))
#' # One runoff file is about 0.5 GB, so this takes a few minutes
#' ws <- npswavacore::get_wsd(park = "ZION")
#' nc <- npswavacore::wava_data_path("all", "loca2_mwbm",
#'                                "runoff_MWBM_LOCA2.ACCESS-CM2.ssp245.r1i1p1f1_2015-2100.nc")
#' out_dir <- file.path(tempdir(), "loca2_runoff")
#' process_loca2_mwbm(nc, ws, var = "runoff",
#'                    historical = c("1979-01-01", "2014-12-31"),
#'                    future = c("2015-01-01", "2100-12-31"),
#'                    output_dir = out_dir, grab_nearest = TRUE)
#' head(arrow::read_parquet(list.files(out_dir, full.names = TRUE)[1]))

process_loca2_mwbm <- function(nc_path,
                               water_supplies,
                               var,
                               historical = c("1980-01-01", "2024-12-31"),
                               future = c("2035-01-01", "2065-12-31"),
                               output_dir,
                               grab_nearest = FALSE,
                               nearest_max_km = 200){

  stopifnot(
    "nearest_max_km must be one positive number" =
      is.numeric(nearest_max_km) && length(nearest_max_km) == 1L &&
        is.finite(nearest_max_km) && nearest_max_km > 0
  )

  supported_vars <- c(
    "aet",
    "deficit",
    "pet",
    "runoff",
    "snow",
    "swe",
    "stor",
    "pr",
    "pr_tavg",
    "tas",
    "tasmin_tavg",
    "tasmax_tavg"
  )

  if (!var %in% supported_vars) {
    stop(
      "Unsupported var '", var, "'. Supported values are: ",
      paste(supported_vars, collapse = ", "), "."
    )
  }
  
  
  message(paste0("Starting vectorized extraction of ", nc_path))
  
  start_clock_time <- Sys.time()
  
  # open using ncdf4 package and ensure it is always closed on exit
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  
  # Some files use different internal NetCDF variable names for the same
  # user-facing variable. Keep the original var label for later output and
  # aggregation decisions; nc_var is what gets passed to ncvar_get.
  nc_var <- var
  
  if(nc_var == "pr" & !grepl("alder", nc_path, ignore.case = TRUE)){
    nc_var <- "pr_tavg"
  }
  
  if(nc_var == "pr_tavg" & grepl("livneh_unsplit_precip", nc_path, ignore.case = TRUE)){
    nc_var <- "PRCP"
  }
  
  if(nc_var == "tasmin_tavg" & grepl("lusu", nc_path, ignore.case = TRUE)){
    nc_var <- "Tmin"
  }
  
  if(nc_var == "tasmax_tavg" & grepl("lusu", nc_path, ignore.case = TRUE)){
    nc_var <- "Tmax"
  }
  
  if(nc_var == "swe"){
    nc_var <- "snow"
  }
  
  # Guard: check that nc_var actually exists in this file before extraction.
  # This makes variable-name mismatches fail early with a useful message.
  available_vars <- names(nc$var)
  if(!nc_var %in% available_vars){
    stop(sprintf(
      "Variable '%s' (mapped from '%s') not found in %s.\nAvailable variables: %s",
      nc_var, var, basename(nc_path), paste(available_vars, collapse = ", ")
    ))
  }
  
  # grab layer of user-selected variable
  v <- nc$var[[nc_var]]
  lon <- ncdf4::ncvar_get(nc, "lon")
  lat <- ncdf4::ncvar_get(nc, "lat")
  
  # If longitudes are 0-360, convert to -180-180
  # Alder training historical data is in 0-360 whereas others are not...
  if (any(lon > 180)) lon <- ifelse(lon > 180, lon - 360, lon)
  
  # transform water_supplies to nc coordinate system if not already
  if (sf::st_crs(water_supplies) != sf::st_crs(terra::rast(nc_path))) {
    water_supplies <- sf::st_transform(water_supplies, sf::st_crs(terra::rast(nc_path)))
  }
  
  # get time units
  # Some Livneh precipitation files use capitalized Time instead of time.
  time_seq <- if(nc_var == "PRCP") {
    ncdf4::ncvar_get(nc, "Time")
  } else {
    ncdf4::ncvar_get(nc, "time")
  }
  
  # Pull the matching time-unit metadata so numeric time values can be
  # converted back to dates below.
  time_units <- if(nc_var == "PRCP") {
    nc$dim[["Time"]]$units
  } else {
    nc$dim[["time"]]$units
  }
  
  start_time_date <- lubridate::ymd(
    sub(" .*", "", sub(".*since\\s+", "", time_units))
  )
  
  # get mwbm units
  mwbm_units <- ncdf4::ncatt_get(nc, nc_var, "units")$value
  
  # Filter time indices based on units
  time_indices_to_keep <- seq_along(time_seq)
  
  if (!is.na(time_units) && grepl("days since", time_units, ignore.case = TRUE)) {
    
    # Convert date boundaries to the same "days since origin" units used by
    # the NetCDF time coordinate.
    origin_date <- as.Date(start_time_date)
    
    # Use the user-specified historical and future windows.
    start_historical <- as.numeric(as.Date(historical[1]) - origin_date)
    end_historical   <- as.numeric(as.Date(historical[2]) - origin_date)
    start_future     <- as.numeric(as.Date(future[1]) - origin_date)
    end_future       <- as.numeric(as.Date(future[2]) - origin_date)
    
    # Grab only layers that fall within the requested reporting periods.
    period_historical <- which(time_seq >= start_historical & time_seq <= end_historical)
    period_future     <- which(time_seq >= start_future & time_seq <= end_future)
    
    time_indices_to_keep <- sort(c(period_historical, period_future))
  }
  
  if (length(time_indices_to_keep) == 0) {
    message("No dates in ", basename(nc_path), " overlap the requested historical/future windows; skipping.")
    return(invisible(NULL))
  }
  
  # update time_seq and nt to reflect filtered data
  time_seq <- time_seq[time_indices_to_keep]
  nt <- length(time_seq)
  
  # convert sf water_supplies to terra vect object and extract coordinates
  v_pts <- terra::vect(water_supplies)
  xy <- terra::crds(v_pts)
  
  # Convert point coordinates to array indices by finding the nearest NetCDF
  # grid-cell center for each water supply.
  lon_index <- purrr::map_int(xy[, 1], ~ which.min(abs(lon - .x)))
  lat_index <- purrr::map_int(xy[, 2], ~ which.min(abs(lat - .x)))
  center_indices <- cbind(lon_index = lon_index, lat_index = lat_index)
  
  # Build a 3x3 neighborhood around each nearest grid cell
  # (-1, 0, +1 in both lon and lat directions).
  offsets <- expand.grid(lon_offset = -1:1, lat_offset = -1:1)
  
  neighbor_list <- vector("list", nrow(center_indices))
  
  for (i in 1:nrow(center_indices)) {
    # Get center longitude and latitude indices for the current point.
    lon_index <- center_indices[i, "lon_index"]
    lat_index <- center_indices[i, "lat_index"]
    
    # Apply offsets while clamping to valid grid boundaries.
    nbr_lons <- pmax(1, pmin(length(lon), lon_index + offsets$lon_offset))
    nbr_lats <- pmax(1, pmin(length(lat), lat_index + offsets$lat_offset))
    
    neighbor_list[[i]] <- data.frame(
      point_id = i,
      lon_index = nbr_lons,
      lat_index = nbr_lats
    )
  }
  
  neighbor_list <- dplyr::bind_rows(neighbor_list) |>
    # Clamping at raster edges can map multiple offsets to the same cell.
    # Keep each available grid cell once rather than weighting edge cells twice.
    dplyr::distinct(point_id, lon_index, lat_index, .keep_all = TRUE) |>
    # Add actual geographic coordinates for each neighborhood cell.
    dplyr::mutate(
      cell_lon = lon[lon_index],
      cell_lat = lat[lat_index]
    )
  
  # Extract dimension names and sizes from the selected NetCDF variable.
  dnames <- purrr::map_chr(v$dim, "name")
  dsizes <- purrr::map_int(v$dim, "len")
  
  # Identify time dimension by looking for "time" in dimension names.
  time_dim <- which(grepl("time", dnames, ignore.case = TRUE))

  # Pre-allocate matrix for extracted values:
  # rows = point-neighborhood cells; columns = selected time steps.
  extract_cells <- function(cell_lookup) {
    n_extractions <- nrow(cell_lookup)
    result_matrix <- matrix(NA_real_, nrow = n_extractions, ncol = nt)

    # Use segment-based reading. Segments are contiguous runs of selected time
    # indices, split into chunks of at most 600 steps to cap memory use.
    if (length(time_indices_to_keep) > 1) {
      time_diffs <- diff(time_indices_to_keep)
      breaks <- which(time_diffs > 1)
    } else {
      breaks <- integer(0)
    }

    segment_starts <- c(1, breaks + 1)
    segment_ends <- c(breaks, length(time_indices_to_keep))
    result_position <- 1

    for (seg in seq_along(segment_starts)) {
      seg_start_pos <- segment_starts[seg]
      seg_end_pos <- segment_ends[seg]
      seg_length <- seg_end_pos - seg_start_pos + 1
      n_subchunks <- ceiling(seg_length / 600)

      for (sub in seq_len(n_subchunks)) {
        sub_start <- seg_start_pos + (sub - 1) * 600
        sub_end <- min(seg_start_pos + sub * 600 - 1, seg_end_pos)
        sub_length <- sub_end - sub_start + 1
        file_indices <- time_indices_to_keep[sub_start:sub_end]
        file_start <- file_indices[1]
        file_end <- file_indices[length(file_indices)]
        file_count <- file_end - file_start + 1
        relative_indices <- file_indices - file_start + 1

        start <- rep(1L, length(dsizes))
        count <- dsizes

        if (length(time_dim) == 1) {
          start[time_dim] <- file_start
          count[time_dim] <- file_count
        }

        arr <- ncdf4::ncvar_get(
          nc,
          nc_var,
          start = start,
          count = count,
          collapse_degen = FALSE
        )

        if (length(time_dim) == 1) {
          if (time_dim == 1) {
            arr <- aperm(arr, c(2, 3, 1))
          } else if (time_dim == 2) {
            arr <- aperm(arr, c(1, 3, 2))
          } else if (time_dim != 3) {
            stop(paste("Unexpected time dimension position:", time_dim))
          }
        } else {
          stop("No time dimension found in NetCDF file")
        }

        for (i in seq_len(n_extractions)) {
          result_matrix[
            i,
            result_position:(result_position + sub_length - 1)
          ] <- arr[
            cell_lookup$lon_index[i],
            cell_lookup$lat_index[i],
            relative_indices
          ]
        }

        result_position <- result_position + sub_length
      }
    }

    result_matrix
  }

  result_matrix <- extract_cells(neighbor_list)

  if (grab_nearest) {
    point_has_data <- rowsum(
      as.integer(rowSums(is.finite(result_matrix)) > 0),
      group = neighbor_list$point_id,
      reorder = FALSE
    )[, 1] > 0
    missing_point_ids <- which(!point_has_data)

    if (length(missing_point_ids) > 0) {
      search_cells <- dplyr::bind_rows(lapply(missing_point_ids, function(point_id) {
        point_lon <- xy[point_id, 1]
        point_lat <- xy[point_id, 2]
        lat_buffer <- nearest_max_km / 111.32
        lon_buffer <- nearest_max_km / (111.32 * max(abs(cos(point_lat * pi / 180)), 1e-6))
        candidate_lons <- which(abs(lon - point_lon) <= lon_buffer)
        candidate_lats <- which(abs(lat - point_lat) <= lat_buffer)

        if (!length(candidate_lons) || !length(candidate_lats)) return(NULL)

        candidates <- expand.grid(
          lon_index = candidate_lons,
          lat_index = candidate_lats
        )
        candidates$distance_m <- geosphere::distHaversine(
          cbind(lon[candidates$lon_index], lat[candidates$lat_index]),
          c(point_lon, point_lat)
        )
        candidates <- candidates[candidates$distance_m <= nearest_max_km * 1000, , drop = FALSE]
        if (!nrow(candidates)) return(NULL)
        candidates$point_id <- point_id
        candidates
      }))

      if (nrow(search_cells) > 0) {
        search_values <- extract_cells(search_cells)
        search_cells$has_data <- rowSums(is.finite(search_values)) > 0
        nearest_centers <- search_cells |>
          dplyr::filter(has_data) |>
          dplyr::group_by(point_id) |>
          dplyr::slice_min(distance_m, n = 1, with_ties = FALSE) |>
          dplyr::ungroup()

        if (nrow(nearest_centers) > 0) {
          fallback_neighbors <- dplyr::bind_rows(lapply(seq_len(nrow(nearest_centers)), function(i) {
            point_id <- nearest_centers$point_id[i]
            center_lon <- nearest_centers$lon_index[i]
            center_lat <- nearest_centers$lat_index[i]
            data.frame(
              point_id = point_id,
              lon_index = pmax(1, pmin(length(lon), center_lon + offsets$lon_offset)),
              lat_index = pmax(1, pmin(length(lat), center_lat + offsets$lat_offset))
            )
          })) |>
            dplyr::distinct(point_id, lon_index, lat_index, .keep_all = TRUE) |>
            dplyr::mutate(
              cell_lon = lon[lon_index],
              cell_lat = lat[lat_index]
            )

          fallback_values <- extract_cells(fallback_neighbors)
          replace_point_ids <- unique(fallback_neighbors$point_id)
          keep_original <- !neighbor_list$point_id %in% replace_point_ids
          neighbor_list <- dplyr::bind_rows(
            neighbor_list[keep_original, , drop = FALSE],
            fallback_neighbors
          )
          result_matrix <- rbind(
            result_matrix[keep_original, , drop = FALSE],
            fallback_values
          )
        }
      }
    }
  }

  n_extractions <- nrow(neighbor_list)

  # ============================================================================
  # Faster results_dt creation
  # rep.int() does not support each =, so use rep() for time_index.
  # result_matrix is unrolled column-major, which matches rep.int ordering
  # (all extraction IDs for time step 1, then time step 2, etc.).
  # ============================================================================
  
  # Convert the neighborhood lookup to data.table and add extraction IDs for
  # joining extracted values back to grid-cell metadata.
  neighbor_dt <- data.table::as.data.table(neighbor_list)
  neighbor_dt[, extraction_id := .I]
  
  # Extract point attributes (non-spatial data) for joining back to results.
  pts_attr <- data.table::as.data.table(sf::st_drop_geometry(water_supplies))
  pts_attr[, point_id := .I]
  
  # Convert matrix results to long format. as.vector(result_matrix) unrolls the
  # matrix by column, so extraction IDs repeat for each time step.
  results_dt <- data.table::data.table(
    extraction_id = rep.int(seq_len(n_extractions), times = nt),
    time_index = rep(seq_len(nt), each = n_extractions),
    value = as.vector(result_matrix)
  )
  
  results_dt[
    neighbor_dt,
    on = "extraction_id",
    `:=`(
      point_id = i.point_id,
      cell_lon = i.cell_lon,
      cell_lat = i.cell_lat
    )
  ]
  
  results_dt[
    pts_attr,
    on = "point_id",
    `:=`(
      wsd_source_id = i.wsd_source_id,
      park_unit = i.park_unit
    )
  ]
  
  results_dt[, `:=`(
    time_units = time_units,
    nc_file = basename(nc_path),
    cell_lon = as.character(cell_lon),
    cell_lat = as.character(cell_lat),
    date = as.Date(start_time_date) + as.numeric(time_seq[time_index]),
    extraction_id = NULL,
    time_index = NULL,
    point_id = NULL
  )]
  
  # Precipitation can be stored in flux units (kg m-2 s-1, equivalent to
  # mm/sec for water). Convert to monthly depth using seconds/day and days/month.
  if(mwbm_units != "mm" & grepl("pr", var)){
    message(mwbm_units)
    
    results_dt[, ym := as.Date(format(date, "%Y-%m-01"))]
    results_dt[, days_in_month := lubridate::days_in_month(ym)]
    
    results_dt <- results_dt |>
      dplyr::mutate(value = value * 86400 * days_in_month)
    
  } else if(mwbm_units == "K"){
    
    # Temperature variables are stored in Kelvin; convert to Celsius.
    message(mwbm_units)
    
    results_dt <- results_dt |>
      dplyr::mutate(value = value - 273.15)
    
  }
  
  # Use original var (not nc_var) for aggregation check so the branch is
  # reliable regardless of what nc_var was renamed to
  # Livneh precipitation is daily, so aggregate daily values to monthly totals.
  # Use original var here because nc_var may have been renamed to PRCP.
  if(var %in% c("pr") & grepl("livneh", nc_path)){
    
    results_dt <- results_dt[
      , date := as.Date(format(date, "%Y-%m-01"))
    ][
      , .(value = if (any(is.finite(value))) {
        sum(value[is.finite(value)])
      } else {
        NA_real_
      }),
      by = .(wsd_source_id, park_unit, time_units, nc_file, cell_lon, cell_lat, date)
    ]
  }
  
  results_dt <- results_dt[
    ,
    .(wsd_source_id, park_unit, time_units, nc_file, cell_lat, cell_lon, date, value)
  ]
  
  # Create output directory if it does not exist.
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Use nc_var for output naming since it reflects the actual variable extracted
  if(nc_var %in% c("Tmin", "Tmax")){
    output_path <- file.path(output_dir, paste0(basename(nc_path), "_", nc_var, ".parquet"))
  } else {
    output_path <- file.path(output_dir, paste0(basename(nc_path), ".parquet"))
  }
  
  arrow::write_parquet(results_dt, output_path)
  
  end_clock_time <- Sys.time()
  
  message(sprintf(
    "Extraction completed in %.2f seconds",
    as.numeric(difftime(end_clock_time, start_clock_time, units = "secs"))
  ))
  
  # nc_close() is now handled by on.exit() at the top of the function
}
