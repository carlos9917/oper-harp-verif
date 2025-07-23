# Some functions used to find extrema and find locations
# of the extrema
# 

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
         latlon <- get_lat_lon_general(arr,j,i)
         #cat("Max here: ",latlon, " with value ",current_val,"\n")
         cat(latlon, " ",current_val,"\n")
          if (current_val > 30) {
          
          browser()
          }
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

# this function works only with DMI data
# projected in the stereographic grid
# of DMIs radar domain
get_lat_lon_dmi <- function(fc_field, i, j) {
  domain <- attr(fc_field, "domain")
  proj_info <- domain$projection
  dx <- domain$dx
  dy <- domain$dy
  x0 <- proj_info$x_0
  y0 <- proj_info$y_0
  
  x <- x0 + (j - 1) * dx
  y <- y0 + (i - 1) * dy
  
  library(sf)
  
  coords <- data.frame(x = x, y = y)
  
  proj4string <- paste0(
    "+proj=stere +lat_0=", proj_info$lat_0,
    " +lon_0=", proj_info$lon_0,
    " +lat_ts=", proj_info$lat_ts,
    " +x_0=", proj_info$x_0,
    " +y_0=", proj_info$y_0,
    " +ellps=", proj_info$ellps
  )
  
  sf_point <- st_as_sf(coords, coords = c("x", "y"), crs = proj4string)
  sf_point_ll <- st_transform(sf_point, crs = 4326)
  
  latlon <- st_coordinates(sf_point_ll)
  
  return(c(lat = latlon[2], lon = latlon[1]))
}



get_lat_lon_general <- function(fc_field, i, j) {
  domain <- attr(fc_field, "domain")
  proj_info <- domain$projection
  
  # Check projection type
  if (proj_info$proj == "latlong") {
    # For lat/lon grids, indices directly correspond to lat/lon
    # Need to calculate based on grid spacing and domain bounds
    
    # Get domain bounds
    SW <- domain$SW  # [lon, lat] of southwest corner
    NE <- domain$NE  # [lon, lat] of northeast corner
    nx <- domain$nx
    ny <- domain$ny
    
    # Calculate grid spacing
    dlon <- (NE[1] - SW[1]) / (nx - 1)
    dlat <- (NE[2] - SW[2]) / (ny - 1)
    
    # Calculate lat/lon for grid indices
    lon <- SW[1] + (j - 1) * dlon
    lat <- SW[2] + (i - 1) * dlat
    
    return(c(lat = lat, lon = lon))
    
  } else {
    # For projected grids (like stereographic)
    library(sf)
    
    dx <- domain$dx
    dy <- domain$dy
    x0 <- proj_info$x_0
    y0 <- proj_info$y_0
    
    # Calculate projected coordinates
    x <- x0 + (j - 1) * dx
    y <- y0 + (i - 1) * dy
    
    # Build projection string dynamically
    proj4string <- paste0("+proj=", proj_info$proj)
    
    # Add parameters that exist
    if (!is.null(proj_info$lat_0)) proj4string <- paste0(proj4string, " +lat_0=", proj_info$lat_0)
    if (!is.null(proj_info$lon_0)) proj4string <- paste0(proj4string, " +lon_0=", proj_info$lon_0)
    if (!is.null(proj_info$lat_ts)) proj4string <- paste0(proj4string, " +lat_ts=", proj_info$lat_ts)
    if (!is.null(proj_info$x_0)) proj4string <- paste0(proj4string, " +x_0=", proj_info$x_0)
    if (!is.null(proj_info$y_0)) proj4string <- paste0(proj4string, " +y_0=", proj_info$y_0)
    if (!is.null(proj_info$ellps)) proj4string <- paste0(proj4string, " +ellps=", proj_info$ellps)
    if (!is.null(proj_info$lon0)) proj4string <- paste0(proj4string, " +lon_0=", proj_info$lon0)
    
    # Create coordinates and transform
    coords <- data.frame(x = x, y = y)
    sf_point <- st_as_sf(coords, coords = c("x", "y"), crs = proj4string)
    sf_point_ll <- st_transform(sf_point, crs = 4326)
    
    latlon <- st_coordinates(sf_point_ll)
    
    return(c(lat = latlon[2], lon = latlon[1]))
  }
}
