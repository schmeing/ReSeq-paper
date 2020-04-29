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

dispersion_fit_sum_truseq_csv <- read_csv(paste0(input_path,"SRR490124_sum_dispersion_fit.csv"), col_types = cols())
dispersion_fit_product_truseq_csv <- read_csv(paste0(input_path,"SRR490124_mult_dispersion_fit.csv"), col_types = cols())

dispersion_fit_sum_nextera_csv <- read_csv(paste0(input_path,"S5L001_sum_dispersion_fit.csv"), col_types = cols())
dispersion_fit_product_nextera_csv <- read_csv(paste0(input_path,"S5L001_mult_dispersion_fit.csv"), col_types = cols())

text_size <- 24
tick_text_size <- 20
dispersion_fit_sum_truseq_csv %>% mutate(Type="Sum", Adapter="TruSeq") %>%
  rbind( dispersion_fit_product_truseq_csv %>% mutate(Type="Product", Adapter="TruSeq") ) %>%
  rbind( dispersion_fit_sum_nextera_csv %>% mutate(Type="Sum", Adapter="Nextera") ) %>%
  rbind( dispersion_fit_product_nextera_csv %>% mutate(Type="Product", Adapter="Nextera") ) %>%
  rename(FragmentLength = insert_length, MeanPrediction=prediction_mean, SampleMean=sample_mean) %>%
  filter(mean_fit == "nbinom") %>%
  select( FragmentLength, MeanPrediction, SampleMean, Type, Adapter) %>%
  mutate(Adapter = ordered(Adapter, levels=c("TruSeq", "Nextera"))) %>%
  ggplot(aes(x=SampleMean, y=MeanPrediction, color=FragmentLength)) +
    geom_point(size=3) +
    geom_abline(slope=1, intercept=0) +
    scale_x_continuous(limits=c(0,0.032), expand = c(0.0005,0.0005)) +
    scale_y_continuous(limits=c(0,0.032), expand = c(0.0005,0.0005)) +
    facet_grid(Adapter ~ Type) +
    labs(x=expression(mean(k[n])),
         y=expression(mean(mu[n]))) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(.01, .99),
          legend.justification = c("left", "top"),
          legend.background = element_rect(fill="transparent", color=NA),
          panel.grid.minor = element_blank())

ggsave(args[2], width=297, height=210, units="mm")
