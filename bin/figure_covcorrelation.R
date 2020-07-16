suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<6){
  stop("Not enough parameters. Need: file1, file2, name1, name2, output, legendPosition(ul,ll,ur,lr)", call.=FALSE)
} else if(length(args)>6){
  stop("Too many parameters. Need: file1, file2, name1, name2, output, legendPosition(ul,ll,ur,lr)", call.=FALSE)
}

calc_correlation <- function(sample1, sample2){
  mean1 <- mean(sample1)
  mean2 <- mean(sample2)
  
  cov <- (sample1-mean1)*(sample2-mean2)
  var1 <- (sample1-mean1)^2
  var2 <- (sample2-mean2)^2
  
  correlation <- sum(cov)/sqrt(sum(var1))/sqrt(sum(var2))
  
  return(correlation)
}

readcovfile <- function(filename){
  return(read_tsv(filename, col_names = FALSE, col_types = cols()))
}

show_correlation <- function(file1, file2, name1, name2, legend_position){
  sample1 <- readcovfile(file1)
  sample2 <- readcovfile(file2)
  
  tmp <- max(max(sample1$X2), max(sample2$X2))
  tmp <- full_join(data.frame(X2=1:tmp), sample1, "X2")
  tmp <- full_join(tmp, sample2, "X2")
  sample1 <- replace_na(tmp$X3.x, 0)
  sample2 <- replace_na(tmp$X3.y, 0)
  
  pear_corr <- calc_correlation(sample1, sample2)
  spear_corr <- cor(sample1, sample2, method="spearman")
  
  text_size = 36
  tick_text_size <- 32
  data.frame(sample1, sample2) %>%
    ggplot(aes(x=sample1, y=sample2)) +
      geom_bin2d(bins=100) +
      ggtitle(sprintf("Pearson = %.2f  Spearman = %.2f", pear_corr, spear_corr)) +
      labs(x = name1, y = name2, fill = "# positions") +
      theme_bw() +
      theme(  plot.title = element_text( size = text_size),
              axis.text.x = element_text( size = tick_text_size),
              axis.text.y = element_text( size = tick_text_size),
              axis.title.x = element_text( size = text_size),
              axis.title.y = element_text( size = text_size),
              strip.text.x = element_text( size = text_size),
              strip.text.y = element_text( size = text_size),
              legend.title=element_text(size=text_size), 
              legend.text=element_text(size=tick_text_size),
              legend.position = c(ifelse(substr(legend_position,2,2) == "r", 0.99, .01), ifelse(substr(legend_position,1,1) == "l", 0.01, .99)),
              legend.justification = c(ifelse(substr(legend_position,2,2) == "r", "right", "left"), ifelse(substr(legend_position,1,1) == "l", "bottom", "top")),
              legend.box.background = element_rect(fill="transparent"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
      guides(fill = guide_colourbar(barheight = 10)) +
      {if(name2 == "Coverage BEAR sim 2")scale_x_continuous(expand=expand_scale(mult=c(0.05,0.1)))}
}

show_correlation(args[1], args[2], args[3], args[4], args[6])

ggsave(args[5], width=297, height=210, units="mm")












