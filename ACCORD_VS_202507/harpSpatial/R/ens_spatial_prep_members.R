#' Prepare the ensemble models to do verification for
#'
#' @param eps_model The name of the EPS model. Maybe expressed as a vector if
#'   more than one EPS model is wanted, or a list for multimodel EPS.
#' @param members_in The input member numbers. If only one EPS is set in
#'   \code{eps_model} then this is a vector. If more than one EPS is set in
#'   \code{eps_model}, or a multimodel EPS is wanted, then this is a list.
#' @param members_out The ouput member numbers. Must be the same form as
#'   members_in. If not passed, members_out is set to the same as members_in.
#' @param lags For reading files from a lagged forecast with members run at
#'   different times, the lag times are set here. The times are expressed as a
#'   character vector, or a named list of character vectors in the case of more
#'   than one model in \code{eps_model}, with a number followed by a letter
#'   giving the units. The avialable units are d, h, m, s for days, hours,
#'   minutes and seconds. The lags argument, if not set to NULL must have
#'   exactly the same dimensions as members_in.
#'
#' @return A nested tibble with columns eps_model, sub_model, member, member_out, lag.
#' 
#' @export 
#' 
#' @examples
#' 
ens_spatial_prep_members <- function(eps_model,
                                     members_in           = seq(0, 9),
                                     members_out          = members_in,
                                     lags                 = NULL
                                     ){

   
  # Sanity checks and organisation of members_in as a list
  
  lags_passed <- !is.null(lags)
  if (!lags_passed) {
    lags <- list()
  }
  
  if (is.list(eps_model)) {
    multimodel <- TRUE
    eps_models <- names(eps_model)
    
    if (!is.list(members_in) |
        !identical(eps_models, names(members_in))) {
      stop(
        paste(
          "For multimodel, members_in must be a list with the",
          "same names as in the eps_model argument.",
          sep = "\n  "
        ),
        call. = FALSE
      )
    }
    
    if (lags_passed && !is.list(lags)) {
      stop(
        paste(
          "If lags are desired, you must treat as multimodel",
          "with the same names as in the eps_model argument",
          sep = "\n"
        ),
        call. = FALSE
      )
    }
    
    for (eps in eps_models) {
      if (!identical(eps_model[[eps]], names(members_in[[eps]]))) {
        stop(
          paste(
            "Model names specified in members_in do not match those in eps_model.",
            paste0(
              "eps_model = ",
              eps,
              ": ",
              paste0(eps_model[[eps]], collapse = ", ")
            ),
            paste0("members_in = ", eps, ": ", paste0(names(
              members_in[[eps]]
            ), collapse = ", ")),
            sep = "\n  "
          ),
          call. = FALSE
        )
      }
      
      if (!identical(names(members_out[[eps]]), names(members_in[[eps]]))) {
        stop(
          paste(
            "Model names specified in members_out do not match those in members_in.",
            paste0("members_in  = ", eps, ": ", paste0(names(
              members_in[[eps]]
            ), collapse = ", ")),
            paste0("members_out = ", eps, ": ", paste0(names(
              members_out[[eps]]
            ), collapse = ", ")),
            sep = "\n "
          ),
          call. = FALSE
        )
      }
      
      if (!lags_passed) {
        lags[[eps]] <- list()
      }
      
      if (lags_passed &&
          !identical(names(lags[[eps]]), names(members_in[[eps]]))) {
        stop(
          paste(
            "Model names specified in lags do not match those in members_in.",
            paste0("members_in = ", eps, ": ", paste0(names(
              members_in[[eps]]
            ), collapse = ", ")),
            paste0("lags       = ", eps, ": ", paste0(names(lags[[eps]]), collapse = ", ")),
            sep = "\n "
          ),
          call. = FALSE
        )
      }
      
      for (sub_eps in names(members_in[[eps]])) {
        if (length(members_out[[eps]][[sub_eps]]) != length(members_in[[eps]][[sub_eps]])) {
          stop(
            paste(
              "Number of members specified in members_out is not the same as in members_in.",
              paste0(
                "members_in = ",
                eps,
                ": ",
                sub_eps,
                ": ",
                length(members_in[[eps]][[sub_eps]]),
                " members"
              ),
              paste0(
                "members_out = ",
                eps,
                ": ",
                sub_eps,
                ": ",
                length(members_out[[eps]][[sub_eps]]),
                " members"
              ),
              sep = "\n "
            ),
            call. = FALSE
          )
        }
        if (lags_passed &&
            length(lags[[eps]][[sub_eps]]) != length(members_in[[eps]][[sub_eps]])) {
          stop(
            paste(
              "Number of members specified in lags is not the same as in members_in.",
              paste0(
                "members_in = ",
                eps,
                ": ",
                sub_eps,
                ": ",
                length(members_in[[eps]][[sub_eps]]),
                " members"
              ),
              paste0(
                "lags       = ",
                eps,
                ": ",
                sub_eps,
                ": ",
                length(lags[[eps]][[sub_eps]]),
                " members"
              ),
              sep = "\n "
            ),
            call. = FALSE
          )
        }
        
        if (!lags_passed) {
          lags[[eps]][[sub_eps]] <-
            rep("0s", length(members_in[[eps]][[sub_eps]]))
        }
        
      }
      
    } # end loop over eps_models
    
  } else {
    # eps_model is not List
    
    multimodel  <- FALSE
    eps_models  <- eps_model
    
    if (length(eps_models) > 1) {
      if (!is.list(members_in) |
          !is.list(members_out) |
          !identical(eps_models, names(members_in)) |
          !identical(eps_models, names(members_out))) {
        stop(
          paste(
            "If more than one eps_model is specified, the members must",
            "be passed as a named list with the names as those specified",
            "in eps_model",
            sep = "\n  "
          ),
          call. = FALSE
        )
      }
      
    } else {
      if (!is.list(members_in)) {
        members_temp               <- list()
        members_temp[[eps_models]] <- members_in
        members_in                 <- members_temp
      }
      if (!is.list(members_out)) {
        members_temp               <- list()
        members_temp[[eps_models]] <- members_out
        members_out                <- members_temp
      }
      if (lags_passed && !is.list(lags)) {
        lags_temp              <- list()
        lags_temp[[eps_model]] <- lags
        lags                   <- lags_temp
      }
      
      if (!identical(eps_models, names(members_in)) |
          !identical(eps_models, names(members_out))) {
        stop(
          paste(
            "If specifying members as a named list for a single eps, the",
            "name in the list for members_in must match that specified",
            "in eps_model.",
            sep = "\n  "
          ),
          call. = FALSE
        )
      }
      
      if (lags_passed && !identical(eps_models, names(lags)))  {
        stop(
          paste(
            "If specifying lags as a named list for a single eps, the",
            "name in the list for lags must match that specified",
            "in eps_model.",
            sep = "\n  "
          ),
          call. = FALSE
        )
      }
      
    }
    
    members_in_temp  <- list()
    members_out_temp <- list()
    lags_temp        <- list()
    for (eps in eps_models) {
      if (length(members_out[[eps]]) != length(members_in[[eps]])) {
        stop(
          paste(
            "Number of members specified in members_out is not the same as in members_in.",
            paste0(
              "members_in  = ",
              eps,
              ": ",
              length(members_in[[eps]]),
              " members"
            ),
            paste0(
              "members_out = ",
              eps,
              ": ",
              length(members_out[[eps]]),
              " members"
            ),
            sep = "\n "
          ),
          call. = FALSE
        )
      }
      if (!lags_passed) {
        lags[[eps]] <- rep("0s", length(members_in[[eps]]))
      }
      if (length(lags[[eps]]) != length(members_in[[eps]])) {
        stop(
          paste(
            "Number of members specified in lags is not the same as in members_in.",
            paste0(
              "members_in = ",
              eps,
              ": ",
              length(members_in[[eps]]),
              " members"
            ),
            paste0("lags       = ", eps, ": ", length(lags[[eps]]), " members"),
            sep = "\n "
          ),
          call. = FALSE
        )
      }
      members_in_temp[[eps]]         <- list()
      members_out_temp[[eps]]        <- list()
      lags_temp[[eps]]               <- list()
      members_in_temp[[eps]][[eps]]  <- members_in[[eps]]
      members_out_temp[[eps]][[eps]] <- members_out[[eps]]
      lags_temp[[eps]][[eps]]        <- lags[[eps]]
    }
    members_in  <- members_in_temp
    members_out <- members_out_temp
    lags        <- lags_temp
    
  } # end handling of inputs related to multiple or single models
  
  
  
  # if (is.null(stations)) {
  #   warning(
  #     "No stations specified. Default station list used.",
  #     call. = FALSE, immediate. = TRUE
  #   )
  #   stations <- get("station_list")
  # }
  
  # initialise interpolation weights
  # if no clim file given, use something from data_files:
  # find first existing file (if none: give an error) and  use that to get domain
  # TODO: maybe for GRIB, we would want to pass a FA climfile for initialisation?
  #       so should we use the same file_format?
  
  # if (!is.null(clim_file)) {
  #   message("Initialising interpolation.")
  #   init <- initialise_interpolation(
  #     file_format = clim_format,
  #     clim_file   = clim_file,
  #     correct_t2m = correct_t2m,
  #     method      = interpolation_method,
  #     use_mask    = use_mask,
  #     stations    = stations
  #   )
  #   init$is_ensemble = TRUE
  # } else {
  #   # just leave it uninitialised for now
  #   init <- list(stations = stations, is_ensemble = TRUE)
  # }
  
  
  ########################### THE ACTUAL WORK STARTS HERE! ##################################
  
  # Convert members_in to a tibble for easier manipulation
  
  members_res <- tibble::tibble(eps_model = names(members_in)) %>%
    dplyr::mutate(sub_model = purrr::map(members_in, names)) %>%
    dplyr::mutate(
      member      = purrr::modify_depth(members_in, 2, `[`),
      members_out = purrr::modify_depth(members_out, 2, `[`),
      lag         = purrr::modify_depth(lags, 2, `[`)
    )
 
  # TODO: correct this if (tidyr_new_interface()) {
  if (TRUE) {
    members_res <- tidyr::unnest(members_res,-tidyr::one_of("eps_model"))
  } else {
    members_res <- tidyr::unnest(members_res)
  }
  
  return(members_res)
  
}
