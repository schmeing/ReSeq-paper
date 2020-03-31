#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  stop("Working directory must be supplied", call.=FALSE)
} else if(length(args)>1){
  stop("Only a single working directory must be supplied", call.=FALSE)
}

setwd(args[1])

library(multipanelfigure)

figure <- multi_panel_figure(width=c(100,64), height=c(43,43,35), row_spacing=3, column_spacing=3, unit="mm", panel_label_type="lower-alpha")
figure %<>% fill_panel("SRR490124-screen.png", row=1, column=1, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("SRR490124-errors.pdf", row=1, column=2, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("S5L001-screen.png", row=2, column=1, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("S5L001-errors.pdf", row=2, column=2, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("rep-screen.png", row=3, column=1:2, scaling = "fit", label_just = "bottom")

save_multi_panel_figure(figure, "figure_syserror.pdf")