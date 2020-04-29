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

maxlike_fit_csv <- read_csv(paste0(input_path,"SRR490124_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="SRR490124")
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"SRR3191692_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="SRR3191692"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="S5L001"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S1L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="S1L001"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S9L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="S9L001"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"ERR2017816_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="ERR2017816"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"ERR3085830_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="ERR3085830"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"ERR1955542_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="ERR1955542"))

disps <- maxlike_fit_csv %>%
  select(Fit, Data, RefSeqBin, InsertLength, DispersionA, DispersionB, FunctionCalls) %>%
  filter(Fit=="nbinom", FunctionCalls < 10000) %>%
  select(-Fit, -FunctionCalls)

medians <- disps %>%
  group_by(Data) %>%
  summarize(DispersionA=median(DispersionA), DispersionB=median(DispersionB)) %>%
  ungroup()

text_size <- 20
tick_text_size <- 16
disps %>%
  ggplot(aes(x=DispersionA, y=DispersionB, color=Data)) +
    geom_point(size=2) +
    geom_point(data=medians, size=6, shape=1) +
    scale_color_manual(values=c("#D92120","#488BC2","#7FB972","#E6642C","#781C81","#D9AD3C","#BBBBBB","#4065B1")) +
    xlab("Parameter A") +
    ylab("Parameter B") +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(0.99, 0.99),
          legend.justification = c(1, 1))

ggsave(args[2], width=297, height=210, units="mm")