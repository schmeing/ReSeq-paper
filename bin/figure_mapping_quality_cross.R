suppressMessages(library(stringr))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<8){
  stop("Output file, real data, load directory, species, sample, crossspecies, crosssample and at least one simulator must be supplied", call.=FALSE)
} else if(length(args)>12){
  stop("Only output file, real data, load directory, species, sample, crossspecies, crosssample and a maximum of five simulators must be supplied", call.=FALSE)
}

output <- args[1]
real_path <- args[2]
input_path <- args[3]
species <- args[4]
sample <- args[5]
crossspecies <- args[6]
crosssample <- args[7]
sims <- args[8:length(args)]

input_csv <- read_csv(real_path, col_types = cols()) %>% mutate(simulator="real")
for( sim in sims ){
  input_csv <- rbind(input_csv, read_csv(paste0(input_path,"/storage_",species,"_",sample,"_",sim,"_",crossspecies,"_",crosssample,"_eval_mapping-bowtie2-s_mapping_qualities.csv"), col_types = cols()) %>% mutate(simulator=sim))
}

nolegend <- grepl("nolegend", output, fixed=TRUE)
sample <- recode(sample, "SRR490124"="Ec-Hi2000-TruSeq", "SRR3191692"="Ec-Hi2500-TruSeq", "SRR3191692_assembly"="Ec-Hi2500-TruSeq-asm",
                         "S5L001"="Ec-Hi4000-Nextera", "S1L001"="Bc-Hi4000-Nextera", "S9L001"="Rs-Hi4000-Nextera",
                         "ERR2017816"="At-HiX-TruSeq", "ERR3085830"="Mm-HiX-Unknown", "ERR1955542"="Hs-HiX-TruSeq",
                         "PRJEB33197"="Hs-Nova-TruSeq", "DRR058060"="Ec-Mi-TruSeq", "PRJNA562949"="At-BGI")
crosssample <- recode(crosssample, "SRR490124"="Ec-Hi2000-TruSeq", "SRR3191692"="Ec-Hi2500-TruSeq", "SRR3191692_assembly"="Ec-Hi2500-TruSeq-asm",
                 "S5L001"="Ec-Hi4000-Nextera", "S1L001"="Bc-Hi4000-Nextera", "S9L001"="Rs-Hi4000-Nextera",
                 "ERR2017816"="At-HiX-TruSeq", "ERR3085830"="Mm-HiX-Unknown", "ERR1955542"="Hs-HiX-TruSeq",
                 "PRJEB33197"="Hs-Nova-TruSeq", "DRR058060"="Ec-Mi-TruSeq", "PRJNA562949"="At-BGI")

bwt2_mapq_points <- c(0,2,30,42)
text_size <- 28
tick_text_size <- 24

input_csv %>%
  group_by(simulator) %>%
  mutate(total=sum(count)) %>%
  filter(unmapped==0) %>%
  select(-unmapped) %>%
  group_by(simulator) %>%
  arrange(desc(mapq)) %>%
  mutate(count=cumsum(count)) %>%
  ungroup() %>%
  rename(Simulator=simulator) %>%
  mutate(Simulator = factor(Simulator, levels = c("real", sims))) %>%
  ggplot(aes(x=mapq, y=count*100/total, color=Simulator, shape=Simulator)) +
    geom_line(size=4) +
    geom_point(aes(x=if_else(mapq %in% bwt2_mapq_points,mapq,NA_real_)), na.rm=TRUE, size=8, stroke=3) +
    scale_color_manual(values=c("#D92120","#488BC2","#7FB972","#E6642C","#781C81","#D9AD3C")) +
    scale_shape_manual(values=c(0,1,2,6,4,5)) +
    scale_x_reverse() +
    xlab("Mapping quality") +
    ylab("% mapped") +
    ggtitle(paste0(sample, " <- ", crosssample)) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          plot.title = element_text( size= text_size, hjust=0.5),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(0.99, 0.01),
          legend.justification = c(1, 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
          {if(nolegend)theme(legend.position = "None")} +
          guides(shape=FALSE)

ggsave(output, width=297, height=210, units="mm")