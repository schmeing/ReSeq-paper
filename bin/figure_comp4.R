#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args)<2){
  stop("Output and at least one pdf must be supplied", call.=FALSE)
} else if(length(args)>5){
  stop("Only up to 4 pdfs are supported", call.=FALSE)
}

library(multipanelfigure)

figure <- multi_panel_figure(rows=2, width=c(81,81), height=120, row_spacing=3, column_spacing=3, unit="mm", panel_label_type="lower-alpha")
for(i in 2:length(args)){
  figure %<>% fill_panel(args[i], scaling = "fit", label_just = "bottom")
}

save_multi_panel_figure(figure, args[1])