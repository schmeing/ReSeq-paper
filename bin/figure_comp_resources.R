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
  mutate(dataset = ordered(dataset, levels=c("SRR490124","SRR3191692","S5L001","S1L001","S9L001","ERR2017816","ERR3085830","ERR1955542", "PRJEB33197", "DRR058060", "PRJNA562949"))) %>%
  mutate(dataset = reorder(dataset, desc(dataset))) %>%
  mutate(dataset = recode(dataset, "SRR490124"="Ec-Hi2000-TruSeq\n(4.6Mb, 449x)", "SRR3191692"="Ec-Hi2500-TruSeq\n(4.6Mb, 1011x)",
                          "S5L001"="Ec-Hi4000-Nextera\n(4.6Mb, 508x)", "S1L001"="Bc-Hi4000-Nextera\n(5.2Mb, 1528x)", "S9L001"="Rs-Hi4000-Nextera\n(4.6Mb, 2901x)",
                          "ERR2017816"="At-HiX-TruSeq\n(120Mb, 31x)", "ERR3085830"="Mm-HiX-Unknown\n(2.8Gb, 43x)", "ERR1955542"="Hs-HiX-TruSeq\n(3.3Gb, 40x)",
                          "PRJEB33197"="Hs-Nova-TruSeq\n(3.3Gb, 131x)", "DRR058060"="Ec-Mi-TruSeq\n(4.6Mb, 176x)", "PRJNA562949"="At-BGI\n(120Mb, 194x)")) %>%
  mutate(simulator = ordered(simulator, levels=c("ReSeq","ART","pIRS","NEAT","BEAR"))) %>%
  mutate(sim_part = recode(sim_part, "training"="Training", "simulation"="Simulation")) %>%
  mutate(computing = ordered(computing, levels=c("single","mixed","threaded","multi"))) %>%
  mutate(computing = recode(computing, "single"="Single", "mixed"="Mixed", "threaded"="Threaded", "multi"="Multi")) %>%
  rename("Simulator" = "simulator", "Computing"="computing") %>%
  group_by(Simulator, Computing, sim_part) %>% mutate(nans = cumsum(is.na(cpu_time))) %>% ungroup() %>%
  unite(group, Simulator, Computing, nans, remove=FALSE)

text_size <- 28
tick_text_size <- 24
time_ticks <- c(10,20,35,60,2*60,5*60,10*60,20*60,35*60,60*60,2*60*60,4*60*60,7*60*60,12*60*60,24*60*60,2*24*60*60,4*24*60*60,7*24*60*60,14*24*60*60,28*24*60*60,2*30*24*60*60,4*30*24*60*60)
time_labels <- c("10s","20s","35s","1m","2m","5m","10m","20m","35m","1h","2h","4h","7h","12h","1d","2d","4d","1w","2w","4w","2m","4m")

finish_plot <- list(
  scale_color_manual(values=c("#488BC2","#7FB972","#E6642C","#781C81","#D9AD3C")),
  scale_shape_manual(values=c(15,17,16,18), drop = FALSE),
  facet_wrap( . ~ sim_part),
  xlab("Dataset"),
  coord_flip(),
  theme_bw(),
  theme(axis.text.x = element_text( size = tick_text_size, angle = 90, vjust=0.5),
        axis.text.y = element_text( size = tick_text_size, lineheight = 0.7),
        axis.title.x = element_text( size = text_size),
        axis.title.y = element_text( size = text_size),
        strip.text.x = element_text( size = text_size),
        strip.text.y = element_text( size = text_size),
        legend.title=element_text(size=text_size), 
        legend.text=element_text(size=tick_text_size),
        legend.justification = c(1, 1),
        legend.position = c(0.99, 0.99),
        panel.grid.minor = element_blank()))

comp_res %>%
  filter(sim_part == "Training") %>%
  ggplot(aes(x=dataset, y=cpu_time, color=Simulator, group=group, shape=Computing)) +
    stat_summary(fun.y=sum, geom="line", na.rm=TRUE, size=4) +
    geom_point(na.rm=TRUE, size=8) +
    scale_y_log10(breaks=time_ticks, labels=time_labels) +
    ylab("CPU time") +
    finish_plot +
    guides(shape=FALSE)

ggsave(paste0(args[2],"/cpu_training.pdf"), width=297, height=210, units="mm")

comp_res %>%
  filter(sim_part == "Simulation") %>%
  ggplot(aes(x=dataset, y=cpu_time, color=Simulator, group=group, shape=Computing)) +
    stat_summary(fun.y=sum, geom="line", na.rm=TRUE, size=4) +
    geom_point(na.rm=TRUE, size=8) +
    scale_y_log10(breaks=time_ticks, labels=time_labels) +
    ylab("CPU time") +
    finish_plot +
    guides(color=FALSE)

ggsave(paste0(args[2],"/cpu_simulation.pdf"), width=297, height=210, units="mm")

comp_res %>%
  filter(sim_part == "Training") %>%
  ggplot(aes(x=dataset, y=elapsed_time, color=Simulator, group=group, shape=Computing)) +
    stat_summary(fun.y=sum, geom="line", na.rm=TRUE, size=4) +
    geom_point(na.rm=TRUE, size=8) +
    scale_y_log10(breaks=time_ticks, labels=time_labels) +
    ylab("Elapsed time") +
    finish_plot +
    theme(legend.position = "none")

ggsave(paste0(args[2],"/elapsed_training.pdf"), width=297, height=210, units="mm")

comp_res %>%
  filter(sim_part == "Simulation") %>%
  ggplot(aes(x=dataset, y=elapsed_time, color=Simulator, group=group, shape=Computing)) +
    stat_summary(fun.y=sum, geom="line", na.rm=TRUE, size=4) +
    geom_point(na.rm=TRUE, size=8) +
    scale_y_log10(breaks=time_ticks, labels=time_labels) +
    ylab("Elapsed time") +
    finish_plot +
    theme(legend.position = "none")

ggsave(paste0(args[2],"/elapsed_simulation.pdf"), width=297, height=210, units="mm")

memory_ticks = c(1,2,4,10,20,40,100,200,400,1024,2*1024,4*1024,10*1024,20*1024,40*1024,100*1024,200*1024)
memory_labels = c("1MB","2MB","4MB","10MB","20MB","40MB","100MB","200MB","400MB","1GB","2GB","4GB","10GB","20GB","40GB","100GB","200GB")

comp_res %>%
  filter(sim_part == "Training") %>%
  ggplot(aes(x=dataset, y=max_memory, color=Simulator, group=group, shape=Computing)) +
    stat_summary(fun.y=sum, geom="line", na.rm=TRUE, size=4) +
    geom_point(na.rm=TRUE, size=8) +
    scale_y_log10(breaks=memory_ticks, labels=memory_labels) +
    ylab("Max. memory") +
    finish_plot +
    theme(legend.position = "none")

ggsave(paste0(args[2],"/memory_training.pdf"), width=297, height=210, units="mm")

comp_res %>%
  filter(sim_part == "Simulation") %>%
  ggplot(aes(x=dataset, y=max_memory, color=Simulator, group=group, shape=Computing)) +
    stat_summary(fun.y=sum, geom="line", na.rm=TRUE, size=4) +
    geom_point(na.rm=TRUE, size=8) +
    scale_y_log10(breaks=memory_ticks, labels=memory_labels) +
    ylab("Max. memory") +
    finish_plot +
    theme(legend.position = "none")

ggsave(paste0(args[2],"/memory_simulation.pdf"), width=297, height=210, units="mm")
