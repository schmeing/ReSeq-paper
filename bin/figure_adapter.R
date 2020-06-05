suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
library(grid)

args = commandArgs(trailingOnly=TRUE)

if(length(args)<2){
  stop("Input and output file must be supplied", call.=FALSE)
} else if(length(args)>2){
  stop("Only input and output file must be supplied", call.=FALSE)
}

adapter_csv <- read_csv(args[1], col_types = cols())

text_size <- 22
tick_text_size <- 14

accepted_levels <- c("Left2","Left1","Start","Right1","Right2","Right3","Right4","Right5-51","Right52","Right53","Right54","Right55","Right56-58","Right59","Right60","Right61")
level_names <- c("-","A","GATCGGAAGA","G","C","A","C","...","T","G","-","-","...","-","-","-")
names(level_names) <- accepted_levels

p <- adapter_csv %>%
  filter(0 == template_segment & min(tries_left) == tries_left) %>%
  select(position, kmer, counts) %>%
  group_by(position) %>%
  slice(1:6) %>%
  mutate(nuc=row_number()) %>%
  ungroup %>%
  mutate(nuc = as.factor(ifelse("Start"==position, 5, nuc))) %>%
  mutate(nuc = recode(nuc, "1"="A", "2"="C", "3"="G", "4"="T", "5"="N")) %>%
  filter(position %in% accepted_levels) %>%
  #mutate( position = factor(position, levels=c(paste0("Left",100:1), "Start", paste0("Right",1:100) )) ) %>%
  mutate( position = factor(position, levels=accepted_levels) ) %>%
  ggplot(aes(y=reorder(kmer, desc(kmer)), x=counts, color=nuc)) +
    geom_point(size=6) +
    scale_color_manual(values=c('#E6642C','#488BC2','#781C81','#B5BD4C','black')) +
    facet_wrap(position ~ ., scales = "free_y", drop=FALSE, labeller = labeller(position = level_names)) +
    ylab("10-mer") +
    xlab("Counts") +
    ggtitle("TruSeq Adapter, Index 7") +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size, angle = 90, vjust=0.5),
          axis.text.y = element_text( size = tick_text_size, angle = 30),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          plot.title = element_text( size= text_size, hjust=0.5),
          strip.text.x = element_text( size = tick_text_size),
          strip.text.y = element_text( size = tick_text_size),
          legend.position = "none",
          panel.grid.minor = element_blank())

g <- ggplot_gtable(ggplot_build(p))

strips <- which(grepl('strip-', g$layout$name))
pal <- c("black","black","black","black",'#B5BD4C','#781C81',"black","black",'#488BC2','#E6642C','#488BC2',"black","black",'#E6642C',"black",'#781C81')

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- pal[i]
}


pdf(args[2], width=11.6929, height=8.26772)
plot(g)
dev.off()

#ggsave(args[2], width=297, height=210, units="mm")