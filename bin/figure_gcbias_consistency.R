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

maxlike_fit_csv <- read_csv(paste0(input_path,"S5L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Ec-Hi4000-Nextera")
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S1L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Bc-Hi4000-Nextera"))
maxlike_fit_csv <- rbind(maxlike_fit_csv, read_csv(paste0(input_path,"S9L001_sum_maxlike_fit.csv"), col_types = cols()) %>% mutate(Data="Rs-Hi4000-Nextera"))

maxlike_fit <- maxlike_fit_csv %>%
  rename(FragmentLength=InsertLength) %>%
  filter(Fit == "nbinom") %>%
  select(Data, FragmentLength, matches("^GCbias*")) %>%
  gather(GC,Bias,-FragmentLength,-Data) %>%
  filter(Bias!=0) %>%
  mutate(GC=as.integer(substr(GC,7,10))) %>%
  arrange(desc(FragmentLength))

Untermedian <- function(x){
  sort(x)[[length(x)/2]]
}

medians <- maxlike_fit %>%
  group_by(Data, GC) %>%
  summarize(Bias=Untermedian(Bias))

text_size <- 20
tick_text_size <- 16
maxlike_fit %>%
  ggplot(aes(x=GC, y=Bias, color=FragmentLength)) +
    geom_point(na.rm=TRUE) +
    geom_point(data=medians, color="red", shape=19) +
    scale_x_continuous(limits = c(0,100), expand = c(0.005,0.005)) +
    scale_y_continuous(limits = c(0,2), expand = c(0.02,0.02)) +
    facet_grid(Data ~ .) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = tick_text_size),
          strip.text.y = element_text( size = tick_text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          panel.grid.minor = element_blank())

ggsave(args[2], width=297, height=210, units="mm")