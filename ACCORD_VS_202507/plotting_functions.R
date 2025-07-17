

#' Plotting fields with ggplot
#'
#' @param data geofield as tibble
#' @param column which is the column of the tibble containing geolist with geofield
#' @return "ggplot" plots fields
plot_panel_field <- function(data,
			     column,
			     domain      = get_domain(data[[column]]),
			     subdomain   = get_domain(data[[column]]),
			     breaks      = c(0, 0.1, 0.2, 0.5, 1., 5., 10., 15.,
					20., 25., 30., 35., 40., 45., 50., 100),
			     palette     = c("#FFFFFE","#00FE96","#00FEC8","#00FEFE",
					 "#00C8FE","#0096FE","#0032FE","#3200FE",
					 "#6400FE","#9600FE","#C800FE","#FA00FE",
					 "#C800C8", "#960096","#FF0000"),
			     limits      = c(min(breaks), max(breaks)),
			     NAcolour    = "white",
			     with_legend = TRUE,
			     only_legend = FALSE,
			     title       = FALSE,
			     subtitle    = FALSE,
			     plot_name   = FALSE,
			     plot_path   = ""
			     ){
        
	# TODO: Compare with harp's plot_field
	require(ggpubr)
	require(ggplotify)

        message("data:")
	print(data)

	countries <- get_map(dom = domain, poly = FALSE)
	verif_box <- meteogrid:::DomainExtent(subdomain)
	boxdim    <- data.frame(xmin = verif_box$x0,
		   		xmax = verif_box$x1,
		   		ymin = verif_box$y0,
		   		ymax = verif_box$y1)

        plt_units <- data$units

        p <- ggplot() +
        geom_georaster(aes(geofield = !!as.name(column)), data) +
        geom_path(aes(x, y), countries) +
	geom_rect(data=boxdim, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), colour="red", alpha=0) + 
        scale_fill_gradientn(plt_units,
                             limits   = limits,
                             colours  = palette,
                             values   = scales::rescale(breaks, to=c(0, 1)),
                             breaks   = breaks,
			     na.value = NAcolour,
			     guide    = "legend") +
        guides(fill = guide_legend(override.aes = list(size = 0.5)),
                                     reverse = TRUE) +
        coord_equal(expand = FALSE) +
        theme_harp_map() + 
	theme(plot.margin=unit(c(0,0,0,0), "cm"))

        if (!with_legend) {
		p <- p + theme(legend.position = "none")
        }
	if (is.character(title)){
		p <- p + labs(title = title)
	}

	if (is.character(subtitle)){
		p <- p + labs(subtitle = subtitle)
	}

        if (only_legend) {
		p <- as_ggplot(get_legend(p))
        }

	if (is.character(plot_name)){
                if (plot_path == "") {
                        plot_path <- paste0(here::here(), "/PLOTS/")
                }
	        	ggsave(paste0(plot_path, "/", plot_name))
	        	message("Saved plot to: ", paste0(plot_path, "/", plot_name))
	}

	return(as.ggplot(p))

}




plot_fields_ob_fc <- function(verif_results,
		       ob_dttm,
		       fc_dttm,
		       lead_time,
		       plot_path = ""){

        ob_dttm_POSIXct <- as.POSIXct(veri_time, "%Y%m%d%H", tz = "UTC")
        fc_dttm_POSIXct <- as.POSIXct(fc_dttm, "%Y%m%d%H", tz = "UTC")

        #### observations ####
        ob_tbl <- tibble::tibble(
                  valid_dttm   = ob_dttm,
                  parameter    = parameter,
                  # lead_time    = lead_time,
                  fcdate       = ob_dttm,
                  units        = ob_units,
                  !!as.name(observation) := geolist(verif_results$ob_field)
                  )

         ob_title <- paste0(parameter, "   ",
                         observation, "   ",
                         ob_dttm_POSIXct
                     )

         ob_plot_name <- paste0(
                                observation, "_",
                                ob_dttm, "_",
                                parameter,
                                ".png"
                         )

         ob_gg <- plot_panel_field(ob_tbl,
                      observation,
                      title     = ob_title,
                      breaks    = breaks,
                      palette   = palette,
                      plot_path = plot_path,
                      plot_name = ob_plot_name
                      )


         #### forecast ####
         fc_tbl <- tibble::tibble(
                      valid_dttm   = ob_dttm,
                      parameter    = parameter,
                      lead_time    = lead_time,
                      fcdate       = fc_dttm,
                      units        = fc_units,
                      !!as.name(fcst_model) := geolist(verif_results$fc_field)
         )

         fc_title <- paste0(parameter, "   ",
                      fcst_model, "   ",
                      fc_dttm_POSIXct, " + ",
                      lead_time, param$acc_unit

         )



         fc_plot_name <- paste0(
                      fcst_model, "_",
                      fc_dttm, "+", lead_time, "_",
                      parameter,
                      ".png"
         )

         fc_gg <- plot_panel_field(fc_tbl,
                      fcst_model,
                      # subdomain = get_domain(ob_tbl[[observation]]),
                      title     = fc_title,
                      breaks    = breaks,
                      palette   = palette,
                      plot_path = plot_path,
                      plot_name = fc_plot_name
         )
         # library(patchwork)
         # combined_plot <- ob_gg / fc_gg

}



