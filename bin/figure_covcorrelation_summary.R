#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  stop("Working directory must be supplied", call.=FALSE)
} else if(length(args)>1){
  stop("Only a single working directory must be supplied", call.=FALSE)
}

setwd(args[1])

library(multipanelfigure)

figure <- multi_panel_figure(rows=3, width=c(39,39,39,39), height=100, row_spacing=3, column_spacing=3, unit="mm", panel_label_type="lower-alpha")
figure %<>% fill_panel("cor1.pdf", row=1, column=1, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor2.pdf", row=1, column=2, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor3.pdf", row=1, column=3, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor4.pdf", row=1, column=4, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor5.pdf", row=2, column=1, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor6.pdf", row=2, column=2, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor7.pdf", row=2, column=3, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor8.pdf", row=2, column=4, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor9.pdf", row=3, column=1, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor10.pdf", row=3, column=2, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor11.pdf", row=3, column=3, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("cor12.pdf", row=3, column=4, scaling = "fit", label_just = "bottom")

save_multi_panel_figure(figure, "figure_covcorrelation.pdf")