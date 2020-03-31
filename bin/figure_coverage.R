#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  stop("Working directory must be supplied", call.=FALSE)
} else if(length(args)>1){
  stop("Only a single working directory must be supplied", call.=FALSE)
}

setwd(args[1])

library(multipanelfigure)

figure <- multi_panel_figure(rows=2, columns=2, width=170, height=120, row_spacing=3, column_spacing=3, unit="mm", panel_label_type="lower-alpha")
figure %<>% fill_panel("SRR490124-cov.pdf", row=1, column=1, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("SRR490124-dup.pdf", row=1, column=2, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("S5L001-cov.pdf", row=2, column=1, scaling = "fit", label_just = "bottom")
figure %<>% fill_panel("S5L001-dup.pdf", row=2, column=2, scaling = "fit", label_just = "bottom")

save_multi_panel_figure(figure, "figure_coverage.pdf")