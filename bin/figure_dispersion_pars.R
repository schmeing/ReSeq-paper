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

maxlike_fit_csv <- read_csv(paste0(input_path,"SRR490124_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Ec-Hi2000-TruSeq")
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"SRR3191692_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Ec-Hi2500-TruSeq"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Ec-Hi4000-Nextera"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S1L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Bc-Hi4000-Nextera"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S9L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Rs-Hi4000-Nextera"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"ERR2017816_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="At-HiX-TruSeq"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"ERR3085830_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Mm-HiX-Unknown"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"ERR1955542_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Hs-HiX-TruSeq"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"PRJEB33197_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Hs-Nova-TruSeq"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"DRR058060_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="EC-Mi-TruSeq"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"PRJNA562949_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="At-BGI"))

disps <- maxlike_fit_csv %>%
  select(Fit, Data, RefSeqBin, InsertLength, DispersionA, DispersionB, FunctionCalls) %>%
  filter(Fit=="nbinom", FunctionCalls < 10000) %>%
  select(-Fit, -FunctionCalls)

medians <- disps %>%
  group_by(Data) %>%
  summarize(DispersionA=median(DispersionA), DispersionB=median(DispersionB)) %>%
  ungroup()

text_size <- 24
tick_text_size <- 20
disps %>%
  mutate(Data = factor(Data, levels=unique(Data))) %>%
  ggplot(aes(x=DispersionA, y=DispersionB, color=Data)) +
    geom_point(size=4, shape=1, stroke=2) +
    geom_point(data=medians, size=8, shape=4, stroke=3, color="white") +
    geom_point(data=medians, size=8, shape=4, stroke=2) +
    scale_color_manual(values=c("#7FB972","#488BC2","#D92120","#781C81",'#B15928','#B5BD4C',"#4065B1",'#664CFF',"#E6642C","#D9AD3C","#BBBBBB",'black')) +
    xlab(paste("Parameter",intToUtf8(0x1D6FCL))) +
    ylab(paste("Parameter",intToUtf8(0x1D6FDL))) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(0.99, 0.01),
          legend.justification = c(1, 0))

ggsave(args[2], width=297, height=210, units="mm", device=cairo_pdf)