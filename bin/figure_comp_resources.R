suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<2){
  stop("Input file and output directory must be supplied", call.=FALSE)
} else if(length(args)>2){
  stop("Only input file and output directory must be supplied", call.=FALSE)
}

comp_res_csv <- read_csv(args[1], col_types = cols())

comp_res <- comp_res_csv %>%
  mutate(dataset = ordered(dataset, levels=c("SRR490124","SRR3191692","S5L001","S1L001","S9L001","ERR2017816","ERR3085830","ERR1955542"))) %>%
  mutate(dataset = recode(dataset, "SRR490124"="SRR490124\n(4.6Mb, 449x)", "SRR3191692"="SRR3191692\n(4.6Mb, 1011x)",
                          "S5L001"="S5L001\n(4.6Mb, 508x)", "S1L001"="S1L001\n(5.2Mb, 1528x)", "S9L001"="S9L001\n(4.6Mb, 2901x)",
                          "ERR2017816"="ERR2017816\n(120Mb, 31x)", "ERR3085830"="ERR3085830\n(2.8Gb, 43x)", "ERR1955542"="ERR1955542\n(3.3Gb, 40x)")) %>%
  mutate(simulator = ordered(simulator, levels=c("ReSeq","ART","pIRS","NEAT"))) %>%
  rename("Simulator" = "simulator")

training_prep <- function(df){
  df %>%
    filter(sim_part == "training") %>%
    mutate(Simulator = recode(Simulator, "pIRS"="pIRS*", "NEAT"="NEAT*"))
}

simulation_prep <- function(df){
  df %>%
    filter(sim_part == "simulation") %>%
    mutate(Simulator = recode(Simulator, "ART"="ART*", "NEAT"="NEAT*"))
}

text_size <- 20
tick_text_size <- 16
time_ticks <- c(10,20,35,60,2*60,5*60,10*60,20*60,35*60,60*60,2*60*60,4*60*60,7*60*60,12*60*60,24*60*60,2*24*60*60,4*24*60*60,7*24*60*60,14*24*60*60,28*24*60*60,8*7*24*60*60,14*7*24*60*60)
time_labels <- c("10s","20s","35s","1m","2m","5m","10m","20m","35m","1h","2h","4h","7h","12h","1d","2d","4d","1w","2w","4w","8w","14w")

finish_plot <- list(
  scale_color_manual(values=c("#D92120","#488BC2","#7FB972","#E6642C")),
  xlab("Dataset"),
  theme_bw(),
  theme(axis.text.x = element_text( size = tick_text_size, angle = 90, vjust=0.5),
        axis.text.y = element_text( size = tick_text_size),
        axis.title.x = element_text( size = text_size),
        axis.title.y = element_text( size = text_size),
        strip.text.x = element_text( size = text_size),
        strip.text.y = element_text( size = text_size),
        legend.title=element_text(size=text_size), 
        legend.text=element_text(size=tick_text_size),
        legend.position = c(0.01, 0.99),
        legend.justification = c(0, 1),
        panel.grid.minor = element_blank()))

comp_res %>%
  training_prep %>%
  ggplot(aes(x=dataset, y=cpu_time, color=Simulator)) +
    geom_point(na.rm=TRUE, size=4) +
    scale_y_log10(breaks=time_ticks, labels=time_labels) +
    ylab("CPU time (training)") +
    finish_plot

ggsave(paste0(args[2],"/cpu_training.pdf"), width=297, height=210, units="mm")

comp_res %>%
  simulation_prep %>%
  ggplot(aes(x=dataset, y=cpu_time, color=Simulator)) +
    geom_point(na.rm=TRUE, size=4) +
    scale_y_log10(breaks=time_ticks, labels=time_labels) +
    ylab("CPU time (simulation)") +
    finish_plot

ggsave(paste0(args[2],"/cpu_simulation.pdf"), width=297, height=210, units="mm")

comp_res %>%
  training_prep %>%
  ggplot(aes(x=dataset, y=elapsed_time, color=Simulator)) +
  geom_point(na.rm=TRUE, size=4) +
  scale_y_log10(breaks=time_ticks, labels=time_labels) +
  ylab("Elapsed time (training)") +
  finish_plot

ggsave(paste0(args[2],"/elapsed_training.pdf"), width=297, height=210, units="mm")

comp_res %>%
  simulation_prep %>%
  ggplot(aes(x=dataset, y=elapsed_time, color=Simulator)) +
  geom_point(na.rm=TRUE, size=4) +
  scale_y_log10(breaks=time_ticks, labels=time_labels) +
  ylab("Elapsed time (simulation)") +
  finish_plot

ggsave(paste0(args[2],"/elapsed_simulation.pdf"), width=297, height=210, units="mm")

memory_ticks = c(1,2,4,10,20,40,100,200,400,1024,2*1024,4*1024,10*1024,20*1024,40*1024,100*1024,200*1024)
memory_labels = c("1MB","2MB","4MB","10MB","20MB","40MB","100MB","200MB","400MB","1GB","2GB","4GB","10GB","20GB","40GB","100GB","200GB")

comp_res %>%
  training_prep %>%
  ggplot(aes(x=dataset, y=max_memory, color=Simulator)) +
  geom_point(na.rm=TRUE, size=4) +
  scale_y_log10(breaks=memory_ticks, labels=memory_labels) +
  ylab("Max. memory (training)") +
  finish_plot

ggsave(paste0(args[2],"/memory_training.pdf"), width=297, height=210, units="mm")

comp_res %>%
  simulation_prep %>%
  ggplot(aes(x=dataset, y=max_memory, color=Simulator)) +
  geom_point(na.rm=TRUE, size=4) +
  scale_y_log10(breaks=memory_ticks, labels=memory_labels) +
  ylab("Max. memory (simulation)") +
  finish_plot

ggsave(paste0(args[2],"/memory_simulation.pdf"), width=297, height=210, units="mm")