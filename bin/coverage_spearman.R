suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<4){
  stop("Not enough parameters.\nPairwise spearman correlation: file1, file2, file3, output.\nSpearman correlation with file 1: file1, file2, file3, file4, output.", call.=FALSE)
} else if(length(args)>5){
  stop("Too many parameters.\nPairwise spearman correlation: file1, file2, file3, output.\nSpearman correlation with file 1: file1, file2, file3, file4, output.", call.=FALSE)
}

readcovfile <- function(filename, i){
  return(read_tsv(filename, col_names = c("sequence","position",paste0("count",i)), col_types = cols()))
}

samples <- readcovfile(args[1], 1)
for( i in 2:(length(args)-1) ){
  samples <- full_join( samples, readcovfile(args[i], i), c("sequence","position") )
}
samples <- samples %>% replace(is.na(.), 0)

if(length(args)==4){
  correlations <- c(cor(samples[["count1"]], samples[["count2"]], method="spearman"),
                    cor(samples[["count1"]], samples[["count3"]], method="spearman"),
                    cor(samples[["count2"]], samples[["count3"]], method="spearman"))
  correlations <- data.frame(sample1=c(args[1],args[1],args[2]), sample2=c(args[2],args[3],args[3]), correlation=correlations)
} else{
  correlations <- c(cor(samples[["count1"]], samples[["count2"]], method="spearman"),
                    cor(samples[["count1"]], samples[["count3"]], method="spearman"),
                    cor(samples[["count1"]], samples[["count4"]], method="spearman"))
  correlations <- data.frame(sample1=c(args[1],args[1],args[1]), sample2=c(args[2],args[3],args[4]), correlation=correlations)
}

write_delim(correlations, args[length(args)], col_names=FALSE)
