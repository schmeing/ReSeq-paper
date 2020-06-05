#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args)<4){
  stop("Output and at three pdfs must be supplied", call.=FALSE)
} else if(length(args)>4){
  stop("Only output and three pdfs are supported", call.=FALSE)
}

library(multipanelfigure)

figure <- multi_panel_figure(rows=2, width=c(45,45,69), height=110, row_spacing=3, column_spacing=3, unit="mm", panel_label_type="lower-alpha")
figure %<>% fill_panel(args[2], row=1:2, column=1:2, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel(args[3], row=1, column=3, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel(args[4], row=2, column=3, scaling = "fit", label_just = "bottom")

save_multi_panel_figure(figure, args[1])