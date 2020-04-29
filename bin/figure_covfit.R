#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  stop("Working directory must be supplied", call.=FALSE)
} else if(length(args)>1){
  stop("Only a single working directory must be supplied", call.=FALSE)
}

setwd(args[1])

library(multipanelfigure)

figure <- multi_panel_figure(width=c(25,25,25,25,25,25), height=c(60,40), row_spacing=3, column_spacing=3, unit="mm", panel_label_type="lower-alpha")
figure %<>% fill_panel("gcbias.pdf", row=1, column=1:3, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("surbias.pdf", row=1, column=4:6, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("samplemean.pdf", row=2, column=1:2, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("dispa.pdf", row=2, column=3:4, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("dispb.pdf", row=2, column=5:6, scaling = "fit", label_just = "bottom")

save_multi_panel_figure(figure, "figure_covfit.pdf")