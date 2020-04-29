suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<1){
  stop("Working directoy parameter missing", call.=FALSE)
} else if(length(args)>1){
  stop("Too many parameters. Only working directory required.", call.=FALSE)
}

setwd(args[1])

# Input is strand of read, type is strand of fragment
sur_csv <- read_csv("surrounding_forward_first.csv", col_types = cols()) %>% mutate(Type="Start-Forward-First")
sur_csv <- rbind(sur_csv, read_csv("surrounding_forward_second.csv", col_types = cols()) %>% mutate(Type="Start-Reverse-Second"))
sur_csv <- rbind(sur_csv, read_csv("surrounding_reverse_first.csv", col_types = cols()) %>% mutate(Type="End-Reverse-First"))
sur_csv <- rbind(sur_csv, read_csv("surrounding_reverse_second.csv", col_types = cols()) %>% mutate(Type="End-Forward-Second"))

text_size = 20
tick_text_size <- 16
sur_csv %>% 
  filter(base != "N") %>%
  group_by(Type, pos) %>%
  mutate( total=sum(count) ) %>%
  ungroup() %>%
  mutate( percent=count*100/total) %>%
  ggplot(aes(x=pos, y=percent, color=Type)) +
    geom_point(size=6, position=position_dodge(width = 0.25)) +
    scale_color_manual(values=c('#E6642C','#488BC2','#781C81','#B5BD4C')) +
    facet_wrap(base ~ .) +
    xlab("Fragment position") + 
    ylab("Frequency [%]") +
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
            legend.position = "top") +
    guides(color=guide_legend(nrow=2))
    
ggsave("surrounding.pdf", width=297, height=210, units="mm")












