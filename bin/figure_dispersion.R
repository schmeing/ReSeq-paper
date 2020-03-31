suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<3){
  stop("Load directory and two output files must be supplied", call.=FALSE)
} else if(length(args)>3){
  stop("Only load directory and two output files must be supplied", call.=FALSE)
}

input_path = args[1]

dispersion_fit_sum_truseq_csv <- read_csv(paste0(input_path,"SRR490124_sum_dispersion_fit.csv"), col_types = cols())
dispersion_fit_sum_nextera_csv <- read_csv(paste0(input_path,"S5L001_sum_dispersion_fit.csv"), col_types = cols())

text_size <- 20
tick_text_size <- 16

dispersion_fit_sum_truseq_csv %>%
  filter(mean_fit == "gcspline") %>%
  rename(FragmentLength = insert_length, Mean=prediction_mean, Dispersion=dispersion_prediction) %>%
  select(FragmentLength, Mean, Dispersion) %>%
  ggplot(aes(x=Mean, y=Mean/Dispersion, color=FragmentLength)) +
    geom_point(na.rm=TRUE) +
    geom_smooth(method='lm', na.rm=TRUE) +
    labs(x=expression(mean(mu[n])),
         y=expression(mean(mu[n])/r)) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(.02, .98),
          legend.justification = c("left", "top"),
          legend.box.background = element_rect(fill="transparent"),
          panel.grid.minor = element_blank())


ggsave(args[2], width=297, height=210, units="mm")

dispersion_fit_sum_nextera_csv %>%
  filter(mean_fit == "gcspline") %>%
  rename(FragmentLength = insert_length, Mean=prediction_mean, Dispersion=dispersion_prediction) %>%
  select(FragmentLength, Mean, Dispersion) %>%
  ggplot(aes(x=Mean, y=Mean/Dispersion, color=FragmentLength)) +
    geom_point(na.rm=TRUE) +
    geom_smooth(method='lm', na.rm=TRUE) +
    ylim(0,2) +
    labs(x=expression(mean(mu[n])),
         y=expression(mean(mu[n])/r)) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(.98, .98),
          legend.justification = c("right", "top"),
          legend.box.background = element_rect(fill="transparent"),
          panel.grid.minor = element_blank())

ggsave(args[3], width=297, height=210, units="mm")