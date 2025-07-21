# CORRECTED SLX Implementation following Sass (2021)

library(harp)
library(harpVis)
library(harpIO)
library(harpSpatial)
library(dplyr)
library(hdf5r)
library(ggplot2)
library(gridExtra)
library(viridis)
library(scico)
library(rlang)

#### DATA LOADING SECTION ####
cat("=== STARTING DATA LOADING ===\n")

# File paths and configuration
veri_time <- "202507011600"
veri_time <- "202507121200"
ob_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/INCAPlus_1h/inca/2025/07/01"
ob_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/INCAPlus_1h/inca/2025/07/12"

fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/CLAEF1k/20250701/00"
fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/CLAEF1k/20250712/00"
fc_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/CLAEF1k"

ob_file <- paste0(ob_file_path,"/INCAPlus_1h_RR_ANA_202507011600.nc")
fc_file <- paste0(fc_file_path,"/CLAEF00+0016:00.grb")

ob_file <- paste0(ob_file_path,"/INCAPlus_1h_RR_ANA_202507121200.nc")
fc_file <- paste0(fc_file_path,"/CLAEF00+0012:00.grb2")

fc_file_template <- "{YYYY}{MM}{DD}/00/{det_model}+{LDT4}:00.grb2" #for the new files

fcst_model <- "CLAEF00"

fc_file_format = "grib"
lead_time <- 12
fc_dttm         <- strftime(strptime(veri_time, "%Y%m%d%H") - (lead_time * 3600), "%Y%m%d%H")
fc_file_name <- generate_filenames(file_path     = fc_file_path,
                                   file_date     = fc_dttm,
                                   lead_time     = lead_time,
                                   file_template = fc_file_template,
                                   det_model     = fcst_model)




fc_file_path_dini  <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/NEA/NEA19082700/tp_NEA1908270903.grib"
# Read model precipitation from GRIB2                  202507020200
precip_fc <- harpIO::read_grid(fc_file_path_dini, "Pcp")

obs_file_path <- "/media/cap/extra_work/verification/oper-harp-verif/ACCORD_VS_202507/sample_data/radar_dmi/sqpe/kavrrad_1h/2025/07"
precip_ob <- harpIO::read_grid(
  paste0(obs_file_path,"/202507020200.kavrRAD.01.h5"),"Pcp",
  hdf5_opts = hdf5_opts(data_path = "/pcp/data1/data", odim = FALSE, meta = TRUE)
) 
dom_ob <-  get_domain(precip_ob)

browser()
precip_fc_regrid <- geo_regrid(precip_fc, dom_ob)

obs_field <- as.array(precip_ob)
fc_field <- as.array(precip_fc_regrid) 
fc_field <- as.array(precip_fc)

source("find_local_extrema.R")
source("plot_maxima_utils.R")
tolerance = 0.0
fc_maxima <- find_local_extrema_sass_corrected(precip_fc, "max", tolerance)

browser()
png("fc_maxima.png", width = 800, height = 600, res = 150)
#plot_field_with_extrema_indices(fc_field,fc_maxima, "Fc maxima")
plot_extrema_points_rotated(fc_field,fc_maxima, "Fc maxima")
dev.off()

browser()
library(sf)
maxima_matrix <- matrix(NA, nrow = nrow(fc_field), ncol = ncol(fc_field))
if (nrow(fc_maxima) > 0) {
  rows <- fc_maxima[, "row"]
  cols <- fc_maxima[, "col"]
  values <- fc_maxima[, "value"]

for (i in seq_along(rows)) {
  cat("Assigning value", values[i], "to position [", rows[i], ",", cols[i], "]\n")
  maxima_matrix[rows[i], cols[i]] <- values[i]
}


}

# Write maxima_matrix to CSV file
write.csv(maxima_matrix, file = "maxima_matrix.csv", row.names = FALSE)
dom <- get_domain(fc_field)
#get the coordinates
# Create vectors of x and y coordinates in projection units
x_coords <- seq(from = 0, by = dom$dx, length.out = dom$nx)
y_coords <- seq(from = 0, by = dom$dy, length.out = dom$ny)

#filter out nan indices
valid_indices <- which(fc_maxima[, "col"] >= 1 & fc_maxima[, "col"] <= length(x_coords) &
                       fc_maxima[, "row"] >= 1 & fc_maxima[, "row"] <= length(y_coords))

fc_maxima_valid <- fc_maxima[valid_indices, ]

x_maxima <- x_coords[fc_maxima_valid[, "col"]]
y_maxima <- y_coords[fc_maxima_valid[, "row"]]



#x_maxima <- x_coords[fc_maxima[, "col"]]
#y_maxima <- y_coords[fc_maxima[, "row"]]
library(sf)

# Define the projection string for stereographic (from dom$projection)
proj_stere <- paste0(
  "+proj=stere +lat_0=", dom$projection$lat_0,
  " +lon_0=", dom$projection$lon_0,
  " +k=1 +x_0=", dom$projection$x_0,
  " +y_0=", dom$projection$y_0,
  " +ellps=", dom$projection$ellps,
  " +units=m +no_defs"
)

# Create an sf object with maxima points in projected coordinates
maxima_sf <- st_as_sf(
  data.frame(x = x_maxima, y = y_maxima),
  coords = c("x", "y"),
  crs = proj_stere
)

# Transform to WGS84 lat/lon
maxima_latlon <- st_transform(maxima_sf, crs = 4326)

# Extract lat/lon coordinates
coords <- st_coordinates(maxima_latlon)

browser()


# Load forecast data
cat("Loading forecast data...\n")
# for the old files
#fc_file_opts <- grib_opts(param_find = setNames(list(list(key = 'indicatorOfParameter', value = 61)), "pcp"))
# for the new files
fc_file_opts     <- grib_opts( param_find = setNames(list(use_grib_shortName('tp')), "pcp"))

#precip_fc <- read_grid(file_name=fc_file, parameter="pcp", dttm="20250701", file_format="grib", file_format_opts = fc_file_opts)
#precip_fc <- read_grid(file_name=fc_file, parameter="pcp", dttm=veri_time, file_format="grib", file_format_opts = fc_file_opts)
#ignoring time to check only
#precip_fc <- read_grid(file_name=fc_file, parameter="pcp", file_format="grib", file_format_opts = fc_file_opts)


#parameter <- "accrr1h"

precip_fc <- read_grid(fc_file_name,
                     parameter        = "pcp", #"accrr1h",
                     is_forecast      = TRUE,
                     dttm             = fc_dttm,
                     file_format      = fc_file_format,
                     file_format_opts = fc_file_opts,
                     param_defs       = list(),
                     lead_time        = lead_time)


dom_fc <- get_domain(precip_fc)
cat("Forecast data loaded successfully\n")

# Load observation data
cat("Loading observation data...\n")
parameter_inca <- "acc1h"
ob_file_opts <- netcdf_opts(proj4_var = "lambert_conformal_conic", param_find = list2(!!parameter_inca := "RR"))
precip_ob <- read_grid(file_name=ob_file, parameter="RR", dttm=veri_time, file_format="netcdf", file_format_opts = ob_file_opts)
dom_ob <- get_domain(precip_ob)
cat("Observation data loaded successfully\n")

# Regrid forecast to observation domain
cat("Regridding forecast to observation domain...\n")
precip_fc_regrid <- geo_regrid(precip_fc, dom_ob)
cat("Regridding completed\n")

# Convert to arrays and handle NA values
cat("Converting to arrays and handling NA values...\n")
obs_field <- as.array(precip_ob)
fc_field <- as.array(precip_fc_regrid)

png("obs_field.png", width = 800, height = 600, res = 150)
plot_field(obs_field)
dev.off()
png("fc_field.png", width = 800, height = 600, res = 150)
plot_field(fc_field)
dev.off()

# Handle NA values
#obs_field[is.na(obs_field)] <- 0
#fc_field[is.na(fc_field)] <- 0

cat(sprintf("Grid dimensions: %d x %d\n", nrow(obs_field), ncol(obs_field)))
cat("Data preparation completed\n")

#### SLX IMPLEMENTATION ####

find_local_extrema_sass_corrected <- function(arr, mode="max", tolerance=0.0) {
  #"""
  #Find local extrema exactly as described in Sass (2021)

  #Key correction: For maxima, only consider non-zero values as potential maxima.
  #This prevents every zero-valued dry point from being classified as a maximum.

  #Key points from paper:
  #- Zero-valued dry areas will often exist... multiple points of zero value
  #  will be automatically selected as minima (NOT maxima)
  #- Default tolerance δ = 0 kg/m²
  #- All selected points contribute with equal weight
  #"""
  rows <- nrow(arr)
  cols <- ncol(arr)
  extrema <- matrix(numeric(0), 0, 3, dimnames=list(NULL, c("row", "col", "value")))

  cat(sprintf("  Scanning %d x %d grid for %s extrema...\n", rows, cols, mode))

  # Progress tracking
  total_points <- (rows-2) * (cols-2)
  progress_interval <- max(1, floor(total_points / 20))
  points_processed <- 0

  # Check each point against its 3x3 neighborhood (excluding boundaries)
  for (i in 2:(rows-1)) {
    for (j in 2:(cols-1)) {
      points_processed <- points_processed + 1

      # Show progress every 5%
      if (points_processed %% progress_interval == 0) {
        progress_pct <- round(100 * points_processed / total_points)
        cat(sprintf("    Progress: %d%% (%d/%d points)\r", progress_pct, points_processed, total_points))
        flush.console()
      }

      current_val <- arr[i, j]
      if (is.na(current_val)) next

      # Get 3x3 neighborhood
      neighborhood <- arr[(i-1):(i+1), (j-1):(j+1)]

      is_extremum <- FALSE
      if (mode == "max") {
        # Local maxima should exclude zeros unless they're true peaks
        # Only points with precipitation > tolerance can be maxima
        if (current_val > tolerance) {
          max_neighbor <- max(neighborhood, na.rm = TRUE)
          is_extremum <- (current_val >= max_neighbor - tolerance) && (current_val == max_neighbor)
        }
      } else {
        # Local minima: points that are <= all neighbors within tolerance
        # Paper explicitly states zeros are automatically selected as minima
        min_neighbor <- min(neighborhood, na.rm = TRUE)
        is_extremum <- (current_val <= min_neighbor + tolerance) && (current_val == min_neighbor)
      }

      if (is_extremum) {
        extrema <- rbind(extrema, c(i, j, current_val))
      }
    }
  }

  cat(sprintf("\n  Found %d %s extrema\n", nrow(extrema), mode))
  return(extrema)
}

get_neighbourhood_extreme_sass <- function(arr, i, j, L, mode="max") {
  #"""
  #Get max/min value in (2L+1)×(2L+1) neighbourhood around point (i,j)
  #Following Sass (2021) equations (1a)-(1d)
  #"""
  rows <- nrow(arr)
  cols <- ncol(arr)

  # Define neighbourhood bounds: [i-L, i+L] × [j-L, j+L]
  i_min <- max(1, i - L)
  i_max <- min(rows, i + L)
  j_min <- max(1, j - L)
  j_max <- min(cols, j + L)

  neighborhood <- arr[i_min:i_max, j_min:j_max]

  if (mode == "max") {
    return(max(neighborhood, na.rm = TRUE))
  } else {
    return(min(neighborhood, na.rm = TRUE))
  }
}

score_function_sass <- function(phi, ob, k=0.1, A=4.0) {
  #"""
  #Exact SLX similarity function from Sass (2021) equations (2a)-(2c), (3a)-(3b)
  #"""
  if (is.na(phi) || is.na(ob)) return(NA)

  if (ob > k) {
    if (phi < ob - k) {
      return(phi / max(ob - k, 1e-9))  # Guard against division by zero
    } else if (phi <= ob) {
      return(1.0)
    } else {  # phi > ob
      return(max(1 - (phi - ob) / (A * ob), 0.0))
    }
  } else {  # ob <= k
    if (phi <= k) {
      return(1.0)
    } else {  # phi > k
      return(max(1 - (phi - k) / (A * k), 0.0))
    }
  }
}

calculate_slx_sass_corrected <- function(obs, forecast, neighbourhood_sizes=c(0, 1, 3, 5, 7, 9), 
                                       tolerance=0.0, k=0.1, A=4.0) {
  #"""
  #Calculate SLX scores following Sass (2021) methodology

  #Key correction: Uses corrected extrema detection that doesn't classify
  #all zeros as maxima, which was causing SLX to decrease with neighbourhood size.

  #Parameters match paper specifications:
  #- tolerance: δ parameter (default ≈ 0 kg/m²)
  #- k: dry threshold (default 0.1 kg/m²)
  #- A: penalty parameter (default 4.0)
  #"""

  cat("\n=== FINDING LOCAL EXTREMA  ===\n")

  # Step 1: Find local extrema 
  cat("Finding observed maxima...\n")
  obs_maxima <- find_local_extrema_sass_corrected(obs, "max", tolerance)
  #png("obs_maxima.png", width = 800, height = 600, res = 150)
  #plot_field_with_extrema_indices(obs,obs_maxima, "Obs maxima")
  #dev.off()

  cat("Finding observed minima...\n")
  obs_minima <- find_local_extrema_sass_corrected(obs, "min", tolerance)
  #png("obs_minima.png", width = 800, height = 600, res = 150)
  #plot_field_with_extrema_indices(obs,obs_minima, "Obs minima")
  #dev.off()

  cat("Finding forecast maxima...\n")
  fc_maxima <- find_local_extrema_sass_corrected(forecast, "max", tolerance)
  #png("fcst_maxima.png", width = 800, height = 600, res = 150)
  #plot_field_with_extrema_indices(forecast,fc_maxima, "Fcst maxima")
  #dev.off()

  cat("Finding forecast minima...\n")
  fc_minima <- find_local_extrema_sass_corrected(forecast, "min", tolerance)
  #png("fcst_minima.png", width = 800, height = 600, res = 150)
  #plot_field_with_extrema_indices(forecast,fc_minima, "Fcst minima")
  #dev.off()

  cat("\n Extrema detection results:\n")
  cat(sprintf("Observed maxima: %d (only non-zero precipitation peaks)\n", nrow(obs_maxima)))
  cat(sprintf("Observed minima: %d (includes zeros as per paper)\n", nrow(obs_minima)))
  cat(sprintf("Forecast maxima: %d (only non-zero precipitation peaks)\n", nrow(fc_maxima)))
  cat(sprintf("Forecast minima: %d (includes zeros as per paper)\n", nrow(fc_minima)))

  results <- list()

  cat("\n=== CALCULATING SLX SCORES FOR DIFFERENT NEIGHBORHOOD SIZES ===\n")

  for (L in neighbourhood_sizes) {
    cat(sprintf("\nProcessing neighborhood size L = %d\n", L))

    # Initialize score vectors
    scores_ob_max <- c()
    scores_ob_min <- c()
    scores_fc_max <- c()
    scores_fc_min <- c()

    # Step 2-4: Calculate component scores following equations (4)-(7)

    # SLX_ob_max: Equation (4) - How well forecast captures observed maxima
    if (nrow(obs_maxima) > 0) {
      cat(sprintf("  Computing SLX_ob_max for %d observed maxima...\n", nrow(obs_maxima)))
      for (idx in 1:nrow(obs_maxima)) {
        if (idx %% max(1, floor(nrow(obs_maxima)/10)) == 0) {
          cat(sprintf("    Processing observed maximum %d/%d\r", idx, nrow(obs_maxima)))
          flush.console()
        }
        i <- obs_maxima[idx, "row"]
        j <- obs_maxima[idx, "col"]
        ob_val <- obs_maxima[idx, "value"]
        fc_neighbourhood_max <- get_neighbourhood_extreme_sass(forecast, i, j, L, "max")
        #if (!is.null(fc_neighbourhood_max)) {
        #    if (NROW(fc_neighbourhood_max) > 1) {
        #    print(fc_neighbourhood_max)
        #    png("fc_max_inter.png", width = 800, height = 600, res = 150)
        #    plot_field_with_extrema_indices(forecast, fc_neighbourhood_max, "FC maxima")
        #    dev.off()
        #      }
        #  }
        score <- score_function_sass(fc_neighbourhood_max, ob_val, k, A)
        scores_ob_max <- c(scores_ob_max, score)
      }
      cat("\n")
    }

    # SLX_ob_min: Equation (5) - How well forecast captures observed minima
    if (nrow(obs_minima) > 0) {
      cat(sprintf("  Computing SLX_ob_min for %d observed minima...\n", nrow(obs_minima)))
      for (idx in 1:nrow(obs_minima)) {
        if (idx %% max(1, floor(nrow(obs_minima)/10)) == 0) {
          cat(sprintf("    Processing observed minimum %d/%d\r", idx, nrow(obs_minima)))
          flush.console()
        }
        i <- obs_minima[idx, "row"]
        j <- obs_minima[idx, "col"]
        ob_val <- obs_minima[idx, "value"]
        fc_neighbourhood_min <- get_neighbourhood_extreme_sass(forecast, i, j, L, "min")
        score <- score_function_sass(fc_neighbourhood_min, ob_val, k, A)
        scores_ob_min <- c(scores_ob_min, score)
      }
      cat("\n")
    }

    # SLX_fc_max: Equation (6) - How well observed field captures forecast maxima
    if (nrow(fc_maxima) > 0) {
      cat(sprintf("  Computing SLX_fc_max for %d forecast maxima...\n", nrow(fc_maxima)))
      for (idx in 1:nrow(fc_maxima)) {
        if (idx %% max(1, floor(nrow(fc_maxima)/10)) == 0) {
          cat(sprintf("    Processing forecast maximum %d/%d\r", idx, nrow(fc_maxima)))
          flush.console()
        }
        i <- fc_maxima[idx, "row"]
        j <- fc_maxima[idx, "col"]
        fc_val <- fc_maxima[idx, "value"]
        obs_neighbourhood_max <- get_neighbourhood_extreme_sass(obs, i, j, L, "max")
        score <- score_function_sass(fc_val, obs_neighbourhood_max, k, A)
        scores_fc_max <- c(scores_fc_max, score)
      }
      cat("\n")
    }

    # SLX_fc_min: Equation (7) - How well observed field captures forecast minima
    if (nrow(fc_minima) > 0) {
      cat(sprintf("  Computing SLX_fc_min for %d forecast minima...\n", nrow(fc_minima)))
      for (idx in 1:nrow(fc_minima)) {
        if (idx %% max(1, floor(nrow(fc_minima)/10)) == 0) {
          cat(sprintf("    Processing forecast minimum %d/%d\r", idx, nrow(fc_minima)))
          flush.console()
        }
        i <- fc_minima[idx, "row"]
        j <- fc_minima[idx, "col"]
        fc_val <- fc_minima[idx, "value"]
        obs_neighbourhood_min <- get_neighbourhood_extreme_sass(obs, i, j, L, "min")
        score <- score_function_sass(fc_val, obs_neighbourhood_min, k, A)
        scores_fc_min <- c(scores_fc_min, score)
      }
      cat("\n")
    }

    # Step 5: Calculate component averages
    slx_ob_max <- if (length(scores_ob_max) > 0) mean(scores_ob_max, na.rm=TRUE) else 0.0
    slx_ob_min <- if (length(scores_ob_min) > 0) mean(scores_ob_min, na.rm=TRUE) else 0.0
    slx_fc_max <- if (length(scores_fc_max) > 0) mean(scores_fc_max, na.rm=TRUE) else 0.0
    slx_fc_min <- if (length(scores_fc_min) > 0) mean(scores_fc_min, na.rm=TRUE) else 0.0

    # Overall SLX score: Equation (8)
    slx_total <- 0.25 * (slx_ob_max + slx_ob_min + slx_fc_max + slx_fc_min)

    cat(sprintf("  Completed L=%d: SLX=%.4f\n", L, slx_total))

    results[[as.character(L)]] <- list(
      SLX = slx_total,
      SLX_ob_max = slx_ob_max,
      SLX_ob_min = slx_ob_min,
      SLX_fc_max = slx_fc_max,
      SLX_fc_min = slx_fc_min,
      n_obs_max = nrow(obs_maxima),
      n_obs_min = nrow(obs_minima),
      n_fc_max = nrow(fc_maxima),
      n_fc_min = nrow(fc_minima),
      scores_ob_max = scores_ob_max,
      scores_ob_min = scores_ob_min,
      scores_fc_max = scores_fc_max,
      scores_fc_min = scores_fc_min,
      obs_maxima = obs_maxima,
      obs_minima = obs_minima,
      fc_maxima = fc_maxima,
      fc_minima = fc_minima
    )
  }

  return(results)
}

#### VALIDATION FUNCTION ####
validate_slx_behavior <- function(slx_results) {
  #"""
  #Validate that SLX shows expected behavior (increasing with neighborhood size)
  #"""
  cat("\n=== VALIDATING SLX BEHAVIOR ===\n")

  L_vals <- as.numeric(names(slx_results))
  SLX_vals <- sapply(slx_results, function(x) x$SLX)

  # Check if SLX generally increases with L
  increasing_trend <- all(diff(SLX_vals) >= -0.01)  # Allow small decreases due to noise

  cat(sprintf("SLX values: %s\n", paste(sprintf("%.3f", SLX_vals), collapse=", ")))
  cat(sprintf("Generally increasing trend: %s\n", if(increasing_trend) "✓ YES" else "✗ NO"))

  # Check extrema counts are reasonable
  r0 <- slx_results[["0"]]
  total_extrema <- r0$n_obs_max + r0$n_obs_min + r0$n_fc_max + r0$n_fc_min
  reasonable_count <- total_extrema < 1000  # Should be much less than grid size

  cat(sprintf("Total extrema count: %d\n", total_extrema))
  cat(sprintf("Reasonable extrema count: %s\n", if(reasonable_count) "✓ YES" else "✗ NO"))

  if (increasing_trend && reasonable_count) {
    cat("\n✓ SLX implementation appears CORRECT!\n")
  } else {
    cat("\n✗ SLX implementation may have issues.\n")
  }
}

#### IMPROVED VISUALIZATIONS ####

plot_field_ggplot <- function(field_data, extrema_max, extrema_min, title,
                             max_color="red", min_color="blue") {
  
  # Convert harp field to data frame for ggplot
  field_df <- as.data.frame(field_data, xy = TRUE)
  print(field_df)
  # Get domain info for coordinate conversion
  domain <- get_domain(field_data)
  print(domain)
  
  # Convert extrema to geographic coordinates
  convert_extrema <- function(extrema, domain) {
    if (nrow(extrema) == 0) return(data.frame())
    
    # Simple linear transformation (adjust based on your projection)
    lon <- domain$SW_lon + (extrema[, "col"] - 1) * domain$dx
    lat <- domain$SW_lat + (extrema[, "row"] - 1) * domain$dy
     print(lat)
    data.frame(
      lon = lon,
      lat = lat,
      value = extrema[, "value"],
      type = rep(c("max", "min")[1], nrow(extrema))
    )
  }
  
  max_df <- convert_extrema(extrema_max, domain)
  if (nrow(max_df) > 0) max_df$type <- "max"
  
  min_df <- convert_extrema(extrema_min, domain)
  if (nrow(min_df) > 0) min_df$type <- "min"
  
  extrema_df <- rbind(max_df, min_df)
  
  # Create ggplot
  p <- ggplot(field_df, aes(x = x, y = y)) +
    geom_raster(aes(fill = value)) +
    scale_fill_viridis_c(name = "Precipitation\n(mm/h)", na.value = "gray95") +
    coord_equal() +
    theme_minimal() +
    labs(title = title, x = "Longitude", y = "Latitude")
  print("extrema points")
  print(extrema_df)
  
  # Add extrema points
  if (nrow(extrema_df) > 0) {
    p <- p + geom_point(data = extrema_df, 
                       aes(x = lon, y = lat, color = type, shape = type),
                       size = 3, stroke = 2) +
      scale_color_manual(values = c("max" = max_color, "min" = min_color),
                        name = "Extrema") +
      scale_shape_manual(values = c("max" = 3, "min" = 4),
                        name = "Extrema")
  }
  
  return(p)
}



#############################

plot_field_with_extrema_indices <- function(field_data, extrema, title, 
                                            point_color = "red", 
                                            point_pch = 3, 
                                            point_cex = 1.5, 
                                            point_lwd = 2) {
  # Plot the field with image(), flipping y-axis to match matrix orientation
  image(
    1:ncol(field_data), 1:nrow(field_data), 
    t(field_data)[, nrow(field_data):1], 
    col = viridis::viridis(20), 
    main = title, 
    xlab = "Column", 
    ylab = "Row", 
    axes = TRUE
  )
  box()
  
  # Overlay extrema points, flipping the row index to match the image orientation
  if (NROW(extrema) > 0) {
    points(
      extrema[, "col"], 
      nrow(field_data) - extrema[, "row"] + 1, 
      col = point_color, 
      pch = point_pch, 
      cex = point_cex, 
      lwd = point_lwd
    )
  }
}

#plot_field_with_extrema_indices <- function(field_data, extrema, title, 
#    color="red") {
#  image(t(field_data[nrow(field_data):1, ]), col = viridis::viridis(20), 
#        main = title, xlab = "Column", ylab = "Row", axes = FALSE)
#  axis(1, at = seq(0, 1, length.out = ncol(field_data)), labels = 1:ncol(field_data))
#  axis(2, at = seq(0, 1, length.out = nrow(field_data)), labels = nrow(field_data):1)
#  box()
#  if (nrow(extrema) > 0) {
#    points((extrema[, "col"] - 1) / (ncol(field_data) - 1),
#           1 - (extrema[, "row"] - 1) / (nrow(field_data) - 1),
#           col = color, pch = 3, cex = 1.5, lwd = 3)
#  }
#
#  #legend("topright", legend = c("Local Max", "Local Min"), 
#  #       col = c(max_color, min_color), pch = c(3, 4), pt.cex = 1.5, pt.lwd = 3, bg = "white")
#}


plot_field_with_extrema <- function(field_data, extrema_max, extrema_min, title, 
                                   max_color="red", min_color="blue") {
  # Create base field plot
  #plot_field(field_data, 
  #           main = title,
  #           legend_label = "Precipitation (mm/h)",
  #           col = viridis(20),
  #           na_colour = "gray95")

  # Add extrema points if they exist
  if (nrow(extrema_max) > 0) {
    points(extrema_max[, "col"], extrema_max[, "row"], 
           col = max_color, pch = 3, cex = 1.5, lwd = 3)
  }

  if (nrow(extrema_min) > 0) {
    points(extrema_min[, "col"], extrema_min[, "row"], 
           col = min_color, pch = 4, cex = 1.5, lwd = 3)
  }

  # Add legend for extrema
  legend("topright", 
         legend = c("Local Max", "Local Min"), 
         col = c(max_color, min_color), 
         pch = c(3, 4), 
         pt.cex = 1.5, 
         pt.lwd = 3,
         bg = "white")
}

#### MAIN EXECUTION ####
doms_to_check <-c(0,1,3,5,7,10,15,20)
doms_to_check <-c(0,1) # 0,5,10,15,20

# Calculate SLX scores 
cat("\n=== STARTING SLX CALCULATION ===\n")
cat("Calculating SLX scores from Sass (2021) implementation...\n")
slx_results <- calculate_slx_sass_corrected(obs_field, fc_field,neighbourhood_sizes=doms_to_check,
                                       tolerance=0.01, k=0.1, A=4.0)

# Validate the implementation
#validate_slx_behavior(slx_results)

# Display results table
cat("\n=== GENERATING RESULTS TABLE ===\n")
cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("SLX RESULTS \n")
cat(paste(rep("=", 80), collapse=""), "\n")
cat(sprintf("%-3s %-8s %-8s %-8s %-8s %-8s %-12s\n",
            "L", "SLX", "ob_max", "ob_min", "fc_max", "fc_min", "total_extrema"))
cat(paste(rep("-", 80), collapse=""), "\n")

for (L in names(slx_results)) {
  r <- slx_results[[L]]
  total_extrema <- r$n_obs_max + r$n_obs_min + r$n_fc_max + r$n_fc_min
  cat(sprintf("%-3s %-8.4f %-8.4f %-8.4f %-8.4f %-8.4f %-12d\n",
              L, r$SLX, r$SLX_ob_max, r$SLX_ob_min, r$SLX_fc_max, r$SLX_fc_min, total_extrema))
}

cat("\n=== CREATING VISUALIZATIONS ===\n")

# VISUALIZATION 1: Original fields with extrema (L=1)
cat("Creating field plots with extrema...\n")
par(mfrow = c(1, 2), mar = c(4, 4, 3, 6))

#png("fields_with_extrema_L1.png", width = 800, height = 600, res = 150)
#L1_results <- slx_results[["1"]]
#plot_field_with_extrema_indices(obs_field, L1_results$obs_maxima, L1_results$obs_minima, "Observed Precipitation with Extrema")
#dev.off()


#plot_field_ggplot(precip_ob, L1_results$obs_maxima, L1_results$obs_minima, 
#                       "Observed Precipitation with Extrema")

#plot_field_with_extrema(precip_ob, L1_results$obs_maxima, L1_results$obs_minima, 
#                       "Observed Precipitation with Extrema")
#plot_field_with_extrema(precip_fc_regrid, L1_results$fc_maxima, L1_results$fc_minima, 
#                       "Forecast Precipitation with Extrema")
#dev.off()

# VISUALIZATION 2: SLX component evolution
cat("Creating SLX evolution plots...\n")
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

L_vals <- as.numeric(names(slx_results))
SLX_vals <- sapply(slx_results, function(x) x$SLX)
SLX_ob_max_vals <- sapply(slx_results, function(x) x$SLX_ob_max)
SLX_ob_min_vals <- sapply(slx_results, function(x) x$SLX_ob_min)
SLX_fc_max_vals <- sapply(slx_results, function(x) x$SLX_fc_max)
SLX_fc_min_vals <- sapply(slx_results, function(x) x$SLX_fc_min)

# Plot 1: Overall SLX score
png("slx_overall.png", width = 800, height = 600, res = 150)
plot(L_vals, SLX_vals, type = "b", pch = 19, col = "darkgreen", lwd = 3,
     xlab = "Neighborhood size (L)", ylab = "SLX Score",
     main = "Overall SLX Score vs Neighborhood Size",
     ylim = c(0, 1), cex = 1.2)
grid(col = "gray", lty = 2)

# Add performance zones
abline(h = 0.8, col = "green", lty = 2, lwd = 2)
abline(h = 0.6, col = "orange", lty = 2, lwd = 2)
abline(h = 0.4, col = "red", lty = 2, lwd = 2)
text(max(L_vals) * 0.7, 0.9, "Excellent", col = "green", cex = 0.8)
text(max(L_vals) * 0.7, 0.7, "Good", col = "orange", cex = 0.8)
text(max(L_vals) * 0.7, 0.5, "Moderate", col = "red", cex = 0.8)
dev.off()

# Plot 2: Component scores
png("slx_component.png", width = 800, height = 600, res = 150)
plot(L_vals, SLX_ob_max_vals, type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "Neighborhood size (L)", ylab = "Component SLX Score",
     main = "SLX Component Scores",
     ylim = c(0, 1))
lines(L_vals, SLX_ob_min_vals, type = "b", pch = 17, col = "blue", lwd = 2)
lines(L_vals, SLX_fc_max_vals, type = "b", pch = 15, col = "orange", lwd = 2)
lines(L_vals, SLX_fc_min_vals, type = "b", pch = 18, col = "purple", lwd = 2)
legend("topright", legend = c("Obs Max", "Obs Min", "FC Max", "FC Min"),
       col = c("red", "blue", "orange", "purple"),
       pch = c(19, 17, 15, 18), lwd = 2, cex = 0.7)
grid(col = "gray", lty = 2)
dev.off()

# Plot 3: Extrema counts
extrema_counts <- sapply(slx_results, function(x) x$n_obs_max + x$n_obs_min + x$n_fc_max + x$n_fc_min)
png("extrema_counts.png", width = 800, height = 600, res = 150)
plot(L_vals, extrema_counts, type = "b", pch = 19, col = "darkblue", lwd = 2,
     xlab = "Neighborhood size (L)", ylab = "Total Number of Extrema",
     main = "Extrema Count vs Neighborhood Size")
grid(col = "gray", lty = 2)
dev.off()

## Plot 4: Score function visualization
#phi_range <- seq(0, 6, 0.1)
#ob_values <- c(0.05, 1.0, 2.0, 4.0)
#colors <- c("purple", "blue", "green", "red")
#
#plot(phi_range, sapply(phi_range, function(x) score_function_sass(x, ob_values[1])), 
#     type = "l", col = colors[1], lwd = 2,
#     xlab = "Forecast Value (mm)", ylab = "Score",
#     main = "SLX Score Function", ylim = c(0, 1))
#
#for (i in 2:length(ob_values)) {
#  lines(phi_range, sapply(phi_range, function(x) score_function_sass(x, ob_values[i])), 
#        col = colors[i], lwd = 2)
#}
#
#abline(v = 0.1, col = "gray", lty = 2, lwd = 1)
#legend("topright", 
#       legend = c(paste("obs =", ob_values), "k = 0.1"),
#       col = c(colors, "gray"), 
#       lty = c(rep(1, 4), 2), 
#       lwd = 2, cex = 0.7)
#grid(col = "gray", lty = 2)

# Reset plotting parameters
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Key findings:\n")
cat(sprintf("- Peak SLX score: %.3f at L=%d\n", max(SLX_vals), L_vals[which.max(SLX_vals)]))
cat(sprintf("- Score improvement from L=0 to optimal: %.3f\n", max(SLX_vals) - SLX_vals[1]))
