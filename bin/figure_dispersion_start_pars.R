suppressMessages(library(stringr))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<4){
  stop("Load directory and three output files must be supplied", call.=FALSE)
} else if(length(args)>4){
  stop("Only load directory and three output file must be supplied", call.=FALSE)
}

input_path = args[1]

maxlike_fit_csv <- read_csv(paste0(input_path,"S5L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=1.0, b=1.0")
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_b0.5_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=1.0, b=0.5"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_b0.75_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=1.0, b=0.75"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_b2.0_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=1.0, b=2.0"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_b10.0_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=1.0, b=10.0"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_a10.0_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=10.0, b=1.0"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_a0.5_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=0.5, b=1.0"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_a0.2_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=0.2, b=1.0"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_a0.1_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=0.1, b=1.0"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_a0.5_b0.5_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=0.5, b=0.5"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_a0.1_b0.5_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=0.1, b=0.5"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S5L001_a0.314_b1.24_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Start="a=0.314, b=1.24"))

disps <- maxlike_fit_csv %>%
  select(Fit, Start, RefSeqBin, InsertLength, DispersionA, DispersionB, FunctionCalls) %>%
  filter(Fit=="nbinom", FunctionCalls < 10000) %>%
  select(-Fit, -FunctionCalls) %>%
  filter(Start %in% c("a=1.0, b=1.0","a=1.0, b=2.0","a=1.0, b=10.0","a=0.1, b=1.0","a=0.1, b=0.5","a=0.314, b=1.24"))

medians <- disps %>%
  group_by(Start) %>%
  summarize(DispersionA=median(DispersionA), DispersionB=median(DispersionB)) %>%
  ungroup()

text_size <- 24
tick_text_size <- 20
color_pallete <- c("#7FB972","#488BC2","#D92120","#781C81",'#B15928','#B5BD4C',"#4065B1",'#664CFF',"#E6642C","#D9AD3C","#BBBBBB",'black')
disps %>%
  ggplot(aes(x=DispersionA, y=DispersionB, color=Start)) +
    geom_point(size=4, shape=1, stroke=2) +
    geom_point(data=medians, size=8, shape=4, stroke=3, color="white") +
    geom_point(data=medians, size=8, shape=4, stroke=2) +
    scale_color_manual(values=color_pallete[c(1,2,4,9,10,11)], labels=(disps[['Start']] %>% unique %>% sort %>% str_replace("a",intToUtf8(0x1D6FCL)) %>% str_replace("b",intToUtf8(0x1D6FDL))) ) +
    xlab(paste("Parameter",intToUtf8(0x1D6FCL))) +
    ylab(paste("Parameter",intToUtf8(0x1D6FDL))) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(0.99, 0.99),
          legend.justification = c(1, 1))

ggsave(args[2], width=297, height=210, units="mm", device=cairo_pdf)


maxlike_fit_csv %>%
  select(Fit, Start, RefSeqBin, InsertLength, DispersionA, DispersionB, FunctionCalls, LogLikelihood) %>%
  filter(Fit=="nbinom") %>%
  mutate(FunctionCallsConverged = ifelse(FunctionCalls<10000, FunctionCalls, NaN), LogLikelihoodConverged = ifelse(FunctionCalls<10000, LogLikelihood, NaN), FunctionCalls = FunctionCalls %% 10000) %>%
  group_by(Start) %>%
  summarize(FunctionCalls=mean(FunctionCalls), FunctionCallsConverged=mean(FunctionCallsConverged, na.rm = TRUE), LogLikelihood=mean(LogLikelihood), LogLikelihoodConverged=mean(LogLikelihoodConverged, na.rm = TRUE)) %>%
  ungroup() %>%
  rename(All=FunctionCalls,Converged=FunctionCallsConverged) %>%
  gather(Fits, FunctionCalls, All, Converged) %>%
  rename(All=LogLikelihood,Converged=LogLikelihoodConverged) %>%
  gather(Fits2, LogLikelihood, All, Converged) %>%
  filter(Fits == Fits2) %>%
  select(-Fits2) %>%
  filter(Fits == "All") %>%
  mutate(Start=str_replace(Start, "a",intToUtf8(0x1D6FCL)), Start=str_replace(Start, "b",intToUtf8(0x1D6FDL))) %>%
  ggplot(aes(x=LogLikelihood, y=FunctionCalls, color=Start)) +
    geom_point(size=4) +
    geom_text_repel(aes(label=Start), size=8) +
    scale_color_manual(values=color_pallete) +
    xlab("Mean log(likelihood)") +
    ylab("Mean function calls") +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          legend.position = "none")

ggsave(args[3], width=297, height=210, units="mm", device=cairo_pdf)


likes <- maxlike_fit_csv %>%
  select(Fit, Start, RefSeqBin, InsertLength, DispersionA, DispersionB, FunctionCalls, LogLikelihood) %>%
  filter(Fit=="nbinom") %>%
  mutate(FunctionCalls = FunctionCalls %% 10000) %>%
  filter(Start %in% c("a=1.0, b=2.0","a=1.0, b=1.0","a=0.5, b=0.5","a=1.0, b=10.0")) %>%
  mutate(Start=str_replace(Start, "a",intToUtf8(0x1D6FCL)), Start=str_replace(Start, "b",intToUtf8(0x1D6FDL)))

mean_values <- likes %>%
  group_by(Start) %>%
  summarize(FunctionCalls=mean(FunctionCalls), LogLikelihood=mean(LogLikelihood)) %>%
  ungroup()

likes %>%
  ggplot(aes(x=LogLikelihood, y=FunctionCalls, color=Start)) +
    geom_point(size=4, shape=1, stroke=2) +
    geom_point(data=mean_values, size=8, shape=4, stroke=3, color="white") +
    geom_point(data=mean_values, size=8, shape=4, stroke=2) +
    geom_text(data=mean_values, aes(label=Start), vjust=1.6, size=8) +
    scale_color_manual(values=color_pallete[c(5,9,10,11)]) +
    xlab("log(likelihood)") +
    ylab("Function calls") +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          legend.position = "none")

ggsave(args[4], width=297, height=210, units="mm", device=cairo_pdf)