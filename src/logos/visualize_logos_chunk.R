plots_order_comp <- paste0(marb, plots_order)
sel_ids <- grep(marb, list_plots)
lst_marb <- list_plots[sel_ids]
list_sel_paths <- list_plots_comp[sel_ids]
sel_fig1 <- list_sel_paths[grep(plots_order_comp[1], lst_marb)]
sel_fig2 <- list_sel_paths[grep(plots_order_comp[2], lst_marb)]
sel_fig3 <- list_sel_paths[grep(plots_order_comp[3], lst_marb)]
sel_fig4 <- list_sel_paths[grep(plots_order_comp[4], lst_marb)]

i <- i+1

marb <- marb_order[i]