#' @title Extract LOCA2 extreme-precipitation metrics at water supplies
#'
#' @description
#' Reads one NetCDF file of the revised LOCA2 (v20240915) precipitation-extremes
#' collection (Alder 2025, <https://doi.org/10.5066/P13FGPRM>) and writes the
#' yearly metrics at each water supply to a parquet file. Run it once per file;
#' [calc_dry_period_exposure()], [calc_wet_period_exposure()],
#' [calc_dry_period_sensitivity()], [calc_wet_period_sensitivity()] and
#' [calc_flood_exposure()] read the folder.
#'
#' Each value is the mean over the 3 x 3 block of grid cells around the
#' supply's nearest cell (the same block as [process_loca2_mwbm()]), ignoring
#' missing cells. Five variables are extracted:
#' - `dry_spell_max_length`: the year's longest run of dry days (CDD), with a
#'   2.54 mm/day threshold;
#' - `wet_spell_max_length`: the year's longest run of wet days (CWD), with a
#'   2.54 mm/day threshold;
#' - `prmax1day`: the year's highest 1-day precipitation (ETCCDI Rx1day), mm;
#' - `prmax5day`: the year's highest 5-day precipitation (ETCCDI Rx5day), mm;
#' - `pr_days_above_nonzero_99th`: the year's number of days above the model's
#'   1951-1980 wet-day 99th percentile.
#'
#' Variables stored in inches in the NetCDF (`units = "in"`, the two depth
#' metrics) are converted to mm.
#'
#' The model, scenario and run are read from the file name, which must follow
#' `MODEL.SCENARIO.RUN.YEARS.LOCA_16thdeg_pr_metrics.nc`, unless
#' `model_metadata` gives them.
#'
#' The same extraction serves other LOCA2 annual-metric collections through
#' `variables` and `model_metadata`: [process_temp_threshold_metrics()] uses it
#' for the temperature-threshold files.
#'
#' The parquet file is `<nc file name>.parquet` in `output_dir`, with one row
#' per water supply and year and columns `wsd_source_id`, `park_unit`,
#' `year`, `model`, `ssp`, `run`, `nc_file` and one column per variable.
#'
#' @param nc_path String. Path to one precipitation-metrics NetCDF file, from
#'   `all/loca_cmip6_ext_precip_products` in the shared data folder. Stops with an
#'   error if the file doesn't exist, lacks any of the variables, or has a name
#'   that can't be parsed. Required.
#' @param water_supplies sf object. The water supplies, as `POINT` geometries
#'   with `wsd_source_id` and `park_unit` columns, e.g. from
#'   [npswavacore::get_wsd()]. Stops with an error if any point is outside the
#'   LOCA2 grid. Required.
#' @param variables Character vector. The NetCDF variables to extract, each
#'   with lon, lat and time dimensions (extra dimensions must have length 1).
#'   Stops with an error if the file lacks any of them. Default: the five
#'   variables listed above.
#' @param model_metadata List or `NULL`. The model, scenario and run to record.
#'   * A list with one non-empty string each for `model`, `ssp` and `run`:
#'     used as given, for files named differently.
#'   * `NULL`: parsed from the file name.
#'
#'   Default: `NULL`.
#' @param grab_nearest Logical. What to do when none of a supply's 3 x 3
#'   cells has data (e.g. an island beyond the LOCA2 land grid).
#'   * `TRUE`: moves the block to the nearest cell with data within
#'     `nearest_max_km` and uses its 3 x 3 neighborhood instead. Supplies with
#'     no data within that distance keep `NA`.
#'   * `FALSE`: keeps the supply with `NA` values.
#'
#'   Default: `FALSE`.
#' @param nearest_max_km Numeric. The farthest, in km, `grab_nearest` looks for
#'   a cell with data. Ignored when `grab_nearest = FALSE`. Default: `200`, as
#'   in [process_loca2_mwbm()]; it reaches Dry Tortugas, about 180 km from the
#'   nearest LOCA2 cell with data.
#' @param output_dir String. The folder for the parquet file, created if it
#'   doesn't exist. Default:
#'   `rapid/inputs/loca_cmip6_precip_periods_water_supplies` in the shared
#'   data folder, see [npswavacore::wava_data_path()].
#' @param overwrite Logical. What to do when the output file already exists.
#'   * `TRUE`: processes the file again and replaces it.
#'   * `FALSE`: skips the file and returns the existing path.
#'
#'   Default: `FALSE`.
#' @param verbose Logical. Whether to print progress.
#'   * `TRUE`: prints the file being read, each variable, and the run time.
#'   * `FALSE`: runs silently.
#'
#'   Default: `TRUE`.
#'
#' @return The path of the parquet file, invisibly.
#' @export
#' @seealso [calc_dry_period_exposure()], [calc_wet_period_exposure()],
#'   [calc_dry_period_sensitivity()], [calc_wet_period_sensitivity()] and
#'   [calc_flood_exposure()], which read the parquet files;
#'   [process_loca2_mwbm()] for the other LOCA2 variables.
#' @examplesIf interactive() && dir.exists(npswavacore::wava_data_path("all", "loca_cmip6_ext_precip_products"))
#' # One file is about 0.7 GB, so this takes a few minutes
#' ws <- npswavacore::get_wsd(park = "ZION")
#' nc_name <- "ACCESS-CM2.historical.r1i1p1f1.1950-2014.LOCA_16thdeg_pr_metrics.nc"
#' nc <- npswavacore::wava_data_path("all", "loca_cmip6_ext_precip_products", "historical", nc_name)
#' out_file <- process_extreme_precip_data(nc, ws, output_dir = tempdir())
#' head(arrow::read_parquet(out_file))

process_extreme_precip_data <- function(
    nc_path,
    water_supplies,
    variables = c("dry_spell_max_length", "wet_spell_max_length",
                  "prmax1day", "prmax5day",
                  "pr_days_above_nonzero_99th"),
    model_metadata = NULL,
    grab_nearest = FALSE,
    nearest_max_km = 200,
    output_dir = npswavacore::wava_data_path("rapid", "inputs", "loca_cmip6_precip_periods_water_supplies/"),
    overwrite = FALSE,
    verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # How this works, in brief. One NetCDF file holds one climate model run: a
  # value for every year and every grid cell (a square about 6 km across)
  # covering the lower 48 states. For each water supply we find the cell it
  # sits in, average that cell with the 8 cells around it, and save one row
  # per supply per year.
  # ---------------------------------------------------------------------------
  #
  # Step 1. Check the inputs before doing any slow work.
  required_variables <- variables

  stopifnot(
    "variables must be a non-empty character vector" =
      is.character(variables) && length(variables) > 0,
    "grab_nearest must be TRUE or FALSE" =
      isTRUE(grab_nearest) || isFALSE(grab_nearest),
    "nearest_max_km must be one positive number" =
      is.numeric(nearest_max_km) && length(nearest_max_km) == 1L &&
        is.finite(nearest_max_km) && nearest_max_km > 0,
    "nc_path must identify one existing file" = length(nc_path) == 1 &&
      file.exists(nc_path),
    "water_supplies must be an sf object" = inherits(water_supplies, "sf"),
    "water_supplies must include wsd_source_id" =
      "wsd_source_id" %in% names(water_supplies),
    "water_supplies must include park_unit" =
      "park_unit" %in% names(water_supplies)
  )

  if (!all(sf::st_geometry_type(water_supplies) %in% "POINT")) {
    stop("water_supplies must contain POINT geometries.")
  }

  # Step 2. Decide where the result is saved. If that file already exists and
  # we weren't asked to redo it, stop here and return the existing file.
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_path <- file.path(output_dir, paste0(basename(nc_path), ".parquet"))

  if (file.exists(output_path) && !overwrite) {
    if (verbose) message("Output already exists; skipping: ", output_path)
    return(invisible(output_path))
  }

  if (verbose) {
    message("Processing LOCA2 metrics from ", basename(nc_path))
    start_time <- Sys.time()
  }

  # Step 3. Open the file and make sure it holds every variable we were asked
  # for. on.exit() closes the file again however the function ends, even after
  # an error.
  nc <- ncdf4::nc_open(nc_path)
  on.exit(ncdf4::nc_close(nc), add = TRUE)

  missing_variables <- setdiff(required_variables, names(nc$var))

  if (length(missing_variables)) {
    stop("NetCDF is missing required variables: ",
      paste(missing_variables, collapse = ", "))
  }

  # Step 4. Work out which climate model, scenario and run this file is. The
  # precipitation files say so in their file name; the temperature files pass
  # it in through model_metadata.
  if (is.null(model_metadata)) {
    # Parse model metadata from the fixed dataset filename convention:
    # MODEL.SCENARIO.RUN.YEARS.LOCA_16thdeg_pr_metrics.nc
    filename_parts <- strsplit(basename(nc_path), "\\.")[[1]]

    if (length(filename_parts) < 3L) {
      stop("Could not parse model, scenario, and run from: ", basename(nc_path))
    }

    model_id <- filename_parts[1]
    ssp_id <- filename_parts[2]
    run_id <- filename_parts[3]
  } else {
    metadata_ok <- is.list(model_metadata) &&
      all(c("model", "ssp", "run") %in% names(model_metadata)) &&
      all(vapply(
        model_metadata[c("model", "ssp", "run")],
        function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x),
        logical(1)
      ))

    if (!metadata_ok) {
      stop("model_metadata must be a list with one non-empty string each for ",
        "model, ssp and run.")
    }

    model_id <- model_metadata$model
    ssp_id <- model_metadata$ssp
    run_id <- model_metadata$run
  }

  # Step 5. Read the grid: where each cell's centre is (longitude and
  # latitude) and which years the file covers. The first variable serves as a
  # template for how the file is laid out; every variable is checked against
  # it later.
  lon <- ncdf4::ncvar_get(nc, "lon")
  lat <- ncdf4::ncvar_get(nc, "lat")

  template_variable <- nc$var[[required_variables[1]]]
  template_dim_names <- vapply(
    template_variable$dim,
    function(x) x$name,
    character(1)
  )
  template_dim_sizes <- vapply(
    template_variable$dim,
    function(x) x$len,
    integer(1)
  )

  lon_dim <- which(tolower(template_dim_names) %in% c("lon", "longitude"))
  lat_dim <- which(tolower(template_dim_names) %in% c("lat", "latitude"))
  time_dim <- which(grepl("time", template_dim_names, ignore.case = TRUE))

  if (length(lon_dim) != 1L || length(lat_dim) != 1L || length(time_dim) != 1L) {
    stop(
      "Could not identify one longitude, latitude, and time dimension for ",
      required_variables[1], "."
    )
  }

  # Time is stored as a number of days since a start date. Turn it into real
  # dates so each layer can be labelled with its year.
  time_name <- template_dim_names[time_dim]
  time_values <- ncdf4::ncvar_get(nc, time_name)
  time_units <- nc$dim[[time_name]]$units

  if (
    is.null(time_units) ||
      length(time_units) != 1L ||
      is.na(time_units) ||
      !grepl("^days since\\s+", time_units, ignore.case = TRUE)
  ) {
    stop("Expected the time coordinate to use 'days since ...' units.")
  }

  time_origin_text <- trimws(
    sub("^.*since\\s+", "", time_units, ignore.case = TRUE)
  )
  time_origin <- as.Date(sub("[ T].*$", "", time_origin_text))
  layer_dates <- time_origin + as.numeric(time_values)

  if (
    length(layer_dates) != template_dim_sizes[time_dim] ||
      anyNA(layer_dates)
  ) {
    stop("Could not recover one annual date for every NetCDF time step.")
  }

  # Step 6. Put the supplies on the same coordinates as the grid. These files
  # count longitude from 0 to 360 degrees east rather than -180 to 180, so
  # supplies in the western hemisphere get 360 added.
  supply_points <- sf::st_transform(water_supplies, 4326)
  point_xy <- sf::st_coordinates(supply_points)

  if (any(lon > 180)) {
    point_xy[, 1] <- ifelse(
      point_xy[, 1] < 0,
      point_xy[, 1] + 360,
      point_xy[, 1]
    )
  }

  # Stop with a clear message if any supply lies beyond the outer edges of the
  # grid altogether.
  lon_half_step <- stats::median(abs(diff(lon))) / 2
  lat_half_step <- stats::median(abs(diff(lat))) / 2
  outside_grid <-
    point_xy[, 1] < min(lon) - lon_half_step |
    point_xy[, 1] > max(lon) + lon_half_step |
    point_xy[, 2] < min(lat) - lat_half_step |
    point_xy[, 2] > max(lat) + lat_half_step

  if (any(outside_grid)) {
    outside_ids <- water_supplies$wsd_source_id[outside_grid]
    stop(
      "Water supplies fall outside the LOCA2 grid: ",
      paste(outside_ids, collapse = ", ")
    )
  }

  # Step 7. Find the cell each supply sits in: the one whose centre is
  # closest. If a supply sits exactly on the line between two cells, take the
  # eastern or southern one, so results match the earlier version of this
  # code.
  center_indices <- cbind(
    lon_index = vapply(
      point_xy[, 1],
      function(x) {
        distance <- abs(lon - x)
        candidates <- which(distance == min(distance))
        candidates[which.max(lon[candidates])]
      },
      integer(1)
    ),
    lat_index = vapply(
      point_xy[, 2],
      function(x) {
        distance <- abs(lat - x)
        candidates <- which(distance == min(distance))
        candidates[which.min(lat[candidates])]
      },
      integer(1)
    )
  )

  # Step 8 (optional). Rescue supplies that have no data around them. A few
  # supplies, such as those on Dry Tortugas, an island beyond the edge of the
  # land data, have no data in any of their 9 cells. With grab_nearest = TRUE,
  # we look for the closest cell that does have data, up to nearest_max_km
  # away, and use that spot instead. Which cells have data is the same every
  # year, so checking one year of the first variable is enough.
  if (grab_nearest) {
    mask_count <- rep(1L, length(template_dim_sizes))
    mask_count[c(lon_dim, lat_dim)] <- template_dim_sizes[c(lon_dim, lat_dim)]
    mask_values <- ncdf4::ncvar_get(
      nc, template_variable,
      start = rep(1L, length(template_dim_sizes)),
      count = mask_count,
      collapse_degen = FALSE
    )
    cell_has_data <- apply(
      array(is.finite(mask_values), dim = mask_count),
      c(lon_dim, lat_dim),
      any
    )

    # For each supply: does any of its 9 cells have data?
    block_has_data <- vapply(seq_len(nrow(center_indices)), function(i) {
      lon_block <- pmax(1L, pmin(length(lon), center_indices[i, "lon_index"] + -1:1))
      lat_block <- pmax(1L, pmin(length(lat), center_indices[i, "lat_index"] + -1:1))
      any(cell_has_data[lon_block, lat_block])
    }, logical(1))

    # Distances are measured using ordinary -180 to 180 longitudes.
    lon_180 <- ifelse(lon > 180, lon - 360, lon)

    for (i in which(!block_has_data)) {
      point_lon <- ifelse(point_xy[i, 1] > 180, point_xy[i, 1] - 360, point_xy[i, 1])
      point_lat <- point_xy[i, 2]
      # Only look at cells inside a box reaching nearest_max_km in each
      # direction (a degree of latitude is about 111 km; a degree of longitude
      # gets shorter away from the equator), and keep the ones with data.
      lat_buffer <- nearest_max_km / 111.32
      lon_buffer <- nearest_max_km /
        (111.32 * max(abs(cos(point_lat * pi / 180)), 1e-6))
      candidates <- expand.grid(
        lon_index = which(abs(lon_180 - point_lon) <= lon_buffer),
        lat_index = which(abs(lat - point_lat) <= lat_buffer)
      )
      candidates <- candidates[
        cell_has_data[cbind(candidates$lon_index, candidates$lat_index)], ,
        drop = FALSE
      ]

      if (nrow(candidates)) {
        # Distance over the Earth's surface, in km, then drop anything beyond
        # the limit.
        candidates$distance_km <- geosphere::distHaversine(
          cbind(lon_180[candidates$lon_index], lat[candidates$lat_index]),
          c(point_lon, point_lat)
        ) / 1000
        candidates <- candidates[candidates$distance_km <= nearest_max_km, , drop = FALSE]
      }

      if (!nrow(candidates)) {
        if (verbose) {
          message("  No cell with data within ", nearest_max_km, " km of ",
            water_supplies$wsd_source_id[i], "; left as NA.")
        }
        next
      }

      # Use the closest cell with data as this supply's centre cell.
      nearest <- candidates[which.min(candidates$distance_km), ]
      center_indices[i, "lon_index"] <- nearest$lon_index
      center_indices[i, "lat_index"] <- nearest$lat_index

      if (verbose) {
        message("  ", water_supplies$wsd_source_id[i], ": no data in its 3 x 3 ",
          "block; using the cell ", round(nearest$distance_km, 1), " km away.")
      }
    }
  }

  # Step 9. List the 9 cells (a 3 x 3 block) around each supply's centre cell.
  # Averaging 9 cells evens out the quirks of any single cell;
  # process_loca2_mwbm() does the same. At the edge of the grid the block is
  # trimmed to stay inside it, and unique() removes cells that would otherwise
  # be counted twice.
  offsets <- expand.grid(lon_offset = -1:1, lat_offset = -1:1)
  neighbor_lookup <- data.table::rbindlist(
    lapply(seq_len(nrow(center_indices)), function(point_id) {
      data.table::data.table(
        point_id = point_id,
        lon_index = pmax(
          1L,
          pmin(
            length(lon),
            center_indices[point_id, "lon_index"] + offsets$lon_offset
          )
        ),
        lat_index = pmax(
          1L,
          pmin(
            length(lat),
            center_indices[point_id, "lat_index"] + offsets$lat_offset
          )
        )
      )
    })
  )
  neighbor_lookup <- unique(
    neighbor_lookup,
    by = c("point_id", "lon_index", "lat_index")
  )

  # Step 10. Make the empty results table: one row for every supply in every
  # year.
  n_points <- nrow(water_supplies)
  n_layers <- length(time_values)

  processed <- data.table::data.table(
    point_id = rep(seq_len(n_points), times = n_layers),
    year = rep(as.integer(format(layer_dates, "%Y")), each = n_points)
  )

  # Step 11. Fill in one column per variable.
  for (variable in required_variables) {
    if (verbose) message("  Extracting ", variable)

    # Check this variable is laid out like the template (same grid, same
    # years), so every value lands on the right cell and year.
    metric_variable <- nc$var[[variable]]
    dim_names <- vapply(
      metric_variable$dim,
      function(x) x$name,
      character(1)
    )
    dim_sizes <- vapply(
      metric_variable$dim,
      function(x) x$len,
      integer(1)
    )
    variable_lon_dim <- which(tolower(dim_names) %in% c("lon", "longitude"))
    variable_lat_dim <- which(tolower(dim_names) %in% c("lat", "latitude"))
    variable_time_dim <- which(grepl("time", dim_names, ignore.case = TRUE))
    extra_dims <- setdiff(
      seq_along(dim_names),
      c(variable_lon_dim, variable_lat_dim, variable_time_dim)
    )

    if (
      length(variable_lon_dim) != 1L ||
        length(variable_lat_dim) != 1L ||
        length(variable_time_dim) != 1L ||
        dim_sizes[variable_lon_dim] != length(lon) ||
        dim_sizes[variable_lat_dim] != length(lat) ||
        dim_sizes[variable_time_dim] != n_layers ||
        any(dim_sizes[extra_dims] != 1L)
    ) {
      stop(
        variable,
        " does not have the expected lon, lat, time, and singleton extra dimensions."
      )
    }

    # A holding table for the raw values: one row per supply-and-cell pair,
    # one column per year.
    n_extractions <- nrow(neighbor_lookup)
    neighbor_values <- matrix(
      NA_real_,
      nrow = n_extractions,
      ncol = n_layers
    )

    # Read the file a block of years at a time rather than all at once, to
    # keep memory use down (these yearly files normally fit in one block).
    # Each block is rearranged to longitude x latitude x year, and the values
    # for every cell we need are copied out.
    max_chunk_size <- 600L
    n_chunks <- ceiling(n_layers / max_chunk_size)

    for (chunk_id in seq_len(n_chunks)) {
      chunk_start <- (chunk_id - 1L) * max_chunk_size + 1L
      chunk_end <- min(chunk_id * max_chunk_size, n_layers)
      chunk_length <- chunk_end - chunk_start + 1L

      start <- rep(1L, length(dim_sizes))
      count <- dim_sizes
      start[variable_time_dim] <- chunk_start
      count[variable_time_dim] <- chunk_length

      metric_array <- ncdf4::ncvar_get(
        nc,
        variable,
        start = start,
        count = count,
        collapse_degen = FALSE
      )

      dimension_order <- c(
        variable_lon_dim,
        variable_lat_dim,
        variable_time_dim,
        extra_dims
      )
      metric_array <- aperm(metric_array, dimension_order)
      dim(metric_array) <- c(length(lon), length(lat), chunk_length)

      for (extraction_id in seq_len(n_extractions)) {
        neighbor_values[
          extraction_id,
          chunk_start:chunk_end
        ] <- metric_array[
          neighbor_lookup$lon_index[extraction_id],
          neighbor_lookup$lat_index[extraction_id],
          seq_len(chunk_length)
        ]
      }
    }

    # Average the 9 cells for each supply and year. Empty cells (for example
    # over the ocean) are left out of the average rather than counted as zero.
    # If all 9 are empty, the result is NA.
    neighbor_values[!is.finite(neighbor_values)] <- NA_real_

    neighbor_sums <- rowsum(replace(neighbor_values, is.na(neighbor_values), 0),
      group = neighbor_lookup$point_id,
      reorder = FALSE)

    neighbor_counts <- rowsum(1L * !is.na(neighbor_values),
      group = neighbor_lookup$point_id,
      reorder = FALSE)

    point_means <- neighbor_sums / neighbor_counts
    point_means[neighbor_counts == 0] <- NA_real_

    # Rainfall amounts are stored in inches, but the rest of this project uses
    # millimetres (1 inch = 25.4 mm). Counts of days are left as they are.
    if (identical(metric_variable$units, "in")) {
      point_means <- point_means * 25.4
    }

    data.table::set(processed,
      j = variable,
      value = as.vector(point_means))

    # Free the memory before moving on to the next variable.
    rm(metric_array, neighbor_values)
    invisible(gc())
  }

  # Step 12. Label every row with the supply's ID and park and with the model,
  # scenario and run, put the columns in a fixed order, and save the table as
  # a parquet file.
  point_attributes <- data.table::as.data.table(
    sf::st_drop_geometry(water_supplies)[
      ,
      c("wsd_source_id", "park_unit"),
      drop = FALSE
    ]
  )
  point_attributes[, point_id := .I]

  processed[
    point_attributes,
    on = "point_id",
    `:=`(
      wsd_source_id = i.wsd_source_id,
      park_unit = i.park_unit
    )
  ]

  processed[
    ,
    `:=`(
      nc_file = basename(nc_path),
      model = model_id,
      ssp = ssp_id,
      run = run_id,
      point_id = NULL
    )
  ]

  data.table::setcolorder(processed, c("wsd_source_id", "park_unit", "year", "model", "ssp", "run", "nc_file", required_variables))

  arrow::write_parquet(processed, output_path)

  if (verbose) {
    elapsed <- Sys.time() - start_time
    message(
      "LOCA2 extraction completed in ",
      round(elapsed, 2), " ", attr(elapsed, "units"),
      ": ", output_path
    )
  }

  invisible(output_path)
}
