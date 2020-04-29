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

maxlike_fit_csv <- read_csv(paste0(input_path,"SRR490124_sum_maxlike_fit.csv"), col_types = cols()) %>% rename(FragmentLength=InsertLength)

maxlike_fit <- maxlike_fit_csv %>%
  mutate(RefSeq=as.factor(RefSeqBin)) %>%
  select(Fit, RefSeq, FragmentLength, matches("*GCbias*")) %>%
  gather(GC,Bias,-FragmentLength,-RefSeq,-Fit) %>%
  filter("Raw" == substr(GC,1,3) | Fit == "nbinom") %>%
  mutate(Fit = ifelse("Raw" != substr(GC,1,3), "normalized", Fit), GC=ifelse("Raw" == substr(GC,1,3), substr(GC,4,13), GC)) %>%
  filter(Bias!=0) %>%
  mutate(GC=as.integer(substr(GC,7,10))) %>%
  mutate(Fit = ordered(Fit, levels=c("poisson","gcspline","nbinom","normalized"))) %>%
  arrange(desc(FragmentLength))

Untermedian <- function(x){
  sort(x)[[length(x)/2]]
}

# In the real method the medians are weighted, this is an approximation for the plot here, but with the values all from the same reference sequence this should be a very good one
medians <-  maxlike_fit %>%
  filter(Fit == "normalized") %>%
  group_by(Fit, GC) %>%
  summarize(Bias=Untermedian(Bias))

knots <- maxlike_fit_csv %>%
  select(Fit, RefSeqBin, FragmentLength, matches("GCknots*")) %>%
  filter(Fit=="gcspline", RefSeqBin==0) %>%
  filter(FragmentLength==Untermedian(FragmentLength)) %>%
  select(Fit, FragmentLength, matches("GCknots*")) %>%
  gather(Knot,GC,-FragmentLength,-Fit) %>%
  mutate(Knot=as.integer(substr(Knot,7,8))) %>%
  mutate(Fit = ordered(Fit, levels=c("poisson","gcspline","nbinom","normalized")))

text_size <- 20
tick_text_size <- 16
maxlike_fit %>%
  ggplot(aes(x=GC, y=Bias, color=FragmentLength)) +
    geom_point(na.rm=TRUE) +
    geom_point(data=medians, color="red", shape=19) +
    geom_vline(data=knots, aes(xintercept=GC, color=FragmentLength)) +
    scale_shape_manual(values=c(19, 15, 17, 18, 3, 4)) +
    scale_x_continuous(limits = c(0,100), expand = c(0.005,0.005)) +
    scale_y_continuous(limits = c(0,2), expand = c(0.02,0.02)) +
    facet_grid(Fit ~ .) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          panel.grid.minor = element_blank())

ggsave(args[2], width=297, height=210, units="mm")