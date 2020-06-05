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

maxlike_csv <- read_csv(paste0(input_path,"SRR490124_logit_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Logit", Data=as.factor("Ec-Hi2000-TruSeq"))
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"SRR490124_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Log", Data="Ec-Hi2000-TruSeq") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"SRR3191692_logit_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Logit", Data="Ec-Hi2500-TruSeq") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"SRR3191692_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Log", Data="Ec-Hi2500-TruSeq") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"S5L001_logit_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Logit", Data="Ec-Hi4000-Nextera") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"S5L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Log", Data="Ec-Hi4000-Nextera") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"S1L001_logit_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Logit", Data="Bc-Hi4000-Nextera") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"S1L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Log", Data="Bc-Hi4000-Nextera") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"S9L001_logit_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Logit", Data="Rs-Hi4000-Nextera") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"S9L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Log", Data="Rs-Hi4000-Nextera") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"ERR2017816_logit_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Logit", Data="At-HiX-TruSeq") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"ERR2017816_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Log", Data="At-HiX-TruSeq") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"ERR3085830_logit_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Logit", Data="Mm-HiX-Unknown") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"ERR3085830_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Log", Data="Mm-HiX-Unknown") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"ERR1955542_logit_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Logit", Data="Hs-HiX-TruSeq") )
maxlike_csv <- rbind(maxlike_csv, read_csv(paste0(input_path,"ERR1955542_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Type="Log", Data="Hs-HiX-TruSeq") )

maxlike <- maxlike_csv %>%
  select(Fit, Type, Data, RefSeqBin, InsertLength, FunctionCalls, LogLikelihood) %>%
  filter(Fit=="nbinom") %>%
  select(-Fit) %>%
  mutate(Converged=(FunctionCalls < 10000), FunctionCalls = FunctionCalls %% 10000)

means <- maxlike %>%
  group_by(Type, Data) %>%
  summarize(FunctionCalls = mean(FunctionCalls), LogLikelihood=mean(LogLikelihood), Converged=mean(Converged))

text_size <- 22
tick_text_size <- 16
maxlike %>%
  mutate(Data = ordered(Data, levels=c("At-HiX-TruSeq","Mm-HiX-Unknown","Hs-HiX-TruSeq","Ec-Hi4000-Nextera","Bc-Hi4000-Nextera","Rs-Hi4000-Nextera","Ec-Hi2000-TruSeq","Ec-Hi2500-TruSeq"))) %>%
  ggplot(aes(x=LogLikelihood/1000, y=FunctionCalls, color=Type, shape=Converged)) +
    geom_point(size=3) +
    geom_point(data=means, mapping=aes(x=LogLikelihood/1000, y=FunctionCalls, color=Type), size=6, shape=1) +
    scale_color_manual(values=c("#488BC2","#E6642C")) +
    xlab("log(likelihood)/1000") +
    ylab("Function calls") +
    facet_wrap(Data ~ ., scales = "free_x", as.table=FALSE) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(.99, 1.05),
          legend.justification = c("right", "top"),
          panel.spacing.x = unit(1, "lines"),
          panel.grid.minor = element_blank())

ggsave(args[2], width=297, height=210, units="mm")
