suppressMessages(library(stringr))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<2){
  stop("Load directory and output file must be supplied", call.=FALSE)
} else if(length(args)>2){
  stop("Only load directory and the output file must be supplied", call.=FALSE)
}

input_path = args[1]

#maxlike_fit_csv <- read_csv(paste0(input_path,"ERR2017816_sum_maxlike_fit.csv"))
#maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"SRR490124_sum_maxlike_fit.csv")))
maxlike_fit_csv <- read_csv(paste0(input_path,"S5L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="S5L001")
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S1L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="S1L001"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S9L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="S9L001"))

text_size <- 20
tick_text_size <- 16
maxlike_fit_csv %>%
  mutate(RefSeq=as.factor(RefSeqBin)) %>%
  select(Fit, Data, RefSeq, InsertLength, matches("Surbias*")) %>%
  filter(Fit=="nbinom") %>%
  select(-Fit) %>%
  gather(Surrounding,Bias,-InsertLength,-RefSeq,-Data) %>%
  mutate(Surrounding=substr(Surrounding,8,100), Surrounding=str_pad(Surrounding, max(nchar(Surrounding)), side="left", pad="0")) %>%
  mutate(Position=as.integer(substr(Surrounding,1,2))-10) %>%
  mutate(Base=substr(Surrounding,3,3)) %>%
  filter(-3<=Position,Position<=11) %>%
  mutate( Position = as.factor(Position) ) %>%
  ggplot(aes(x=Position, y=Bias, fill=Data, color=Base)) +
    geom_boxplot( position = position_dodge(width=0.8, preserve="single"), width=3 ) +
    geom_vline(xintercept=seq(0.5, 30.5, 1), color="lightgrey") +
    scale_x_discrete(expand = c(0.03,0.03)) +
    scale_color_manual(values=c('#E6642C','#488BC2','#781C81','#B5BD4C')) +
    scale_fill_manual(values=c("white", "grey", "black")) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_line(colour = "lightgrey"))


ggsave(args[2], width=297, height=210, units="mm")