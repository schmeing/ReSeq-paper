suppressMessages(library(stringr))
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<3){
  stop("Load directory, sample and output prefix must be supplied", call.=FALSE)
} else if(length(args)>3){
  stop("Only load directory, sample and output prefix must be supplied", call.=FALSE)
}

input_path = args[1]
sample = args[2]

input_csv <- read_csv(paste0(input_path,"/storage_ecoli_",sample,"_ReSeq_eval_mapping-bowtie2-lowerror_mapping_correctness.csv"), col_types = cols()) %>% mutate(simulator="ReSeq-bowtie2", mapper="bowtie2")
input_csv <- rbind(input_csv, read_csv(paste0(input_path,"/storage_ecoli_",sample,"_ReSeq_bwa_eval_mapping-bowtie2-lowerror_mapping_correctness.csv"), col_types = cols()) %>% mutate(simulator="ReSeq-bwa", mapper="bowtie2"))
input_csv <- rbind(input_csv, read_csv(paste0(input_path,"/storage_ecoli_",sample,"_ReSeq_eval_mapping-bwa-lowerror_mapping_correctness.csv"), col_types = cols()) %>% mutate(simulator="ReSeq-bowtie2", mapper="bwa"))
input_csv <- rbind(input_csv, read_csv(paste0(input_path,"/storage_ecoli_",sample,"_ReSeq_bwa_eval_mapping-bwa-lowerror_mapping_correctness.csv"), col_types = cols()) %>% mutate(simulator="ReSeq-bwa", mapper="bwa"))

nolegend <- grepl("nolegend", args[3], fixed=TRUE)

text_size <- 24
tick_text_size <- 20
input_csv %>%
  group_by(simulator, mapper, negative) %>%
  mutate(maxTP=sum(count)) %>%
  filter(unmapped==0) %>%
  group_by(simulator, mapper, mapq) %>%
  select(-unmapped, -negative) %>%
  mutate(total=sum(count)) %>%
  ungroup() %>%
  filter(correctness!=0) %>%
  mutate(correctness = recode(as.character(correctness), "1"="overlapping", "2"="correct_start", "3"="correct")) %>%
  spread(correctness, count, fill=0) %>%
  mutate(correct_start=correct+correct_start, overlapping=correct_start+overlapping) %>%
  gather(correctness, TP, correct, correct_start, overlapping) %>%
  mutate(FP = total-TP) %>%
  select(-total) %>%
  group_by(simulator,mapper, correctness) %>%
  arrange(desc(mapq)) %>%
  mutate(TP=cumsum(TP)*100/maxTP, FP=cumsum(FP)*100/maxTP) %>%
  ungroup() %>%
  rename(Simulator=simulator, Mapper=mapper, Correctness=correctness) %>%
  arrange(desc(TP), desc(Correctness), desc(mapq), Mapper) %>%
  ggplot(aes(x=FP, y=TP, color=Mapper, shape=Simulator, fill=Correctness)) +
  geom_line(size=4) +
  geom_point(size=8, stroke=8) +
  scale_shape_manual(values=c(21, 22)) +
  scale_color_manual(values=c("#009E73","#0072B2")) +
  scale_fill_manual(values=c("#000000","#888888","#FFFFFF")) +
  xlab("FP [% of simulated positives]") +
  ylab("TP [% of simulated positives]") +
  theme_bw() +
  theme(axis.text.x = element_text( size = tick_text_size),
        axis.text.y = element_text( size = tick_text_size),
        axis.title.x = element_text( size = text_size),
        axis.title.y = element_text( size = text_size),
        strip.text.x = element_text( size = text_size),
        strip.text.y = element_text( size = text_size),
        legend.title=element_text(size=text_size), 
        legend.text=element_text(size=tick_text_size),
        legend.position = c(0.99, 0.01),
        legend.justification = c(1, 0),
        legend.box = "vertical",
        legend.direction = "vertical",
        legend.key.size = unit(1.5, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = 21), ncol = 3),
         color = guide_legend(ncol = 2),
         shape = guide_legend(ncol = 2)) +
  {if(nolegend)theme(legend.position = "None")}

ggsave(paste0(args[3],"_correctness.pdf"), width=297, height=210, units="mm")

input_real_csv <- read_csv(paste0(input_path,"/storage_ecoli_",sample,"_real_mapping-bowtie2-s_mapping_qualities.csv"), col_types = cols()) %>% mutate(simulator="real", mapper="bowtie2")
input_real_csv <- rbind(input_real_csv, read_csv(paste0(input_path,"/storage_ecoli_",sample,"_real_mapping-bwa-s_mapping_qualities.csv"), col_types = cols()) %>% mutate(simulator="real", mapper="bwa"))

input_real <- input_real_csv %>%
  group_by(simulator, mapper) %>%
  mutate(total=sum(count)) %>%
  filter(unmapped==0) %>%
  select(-unmapped) %>%
  arrange(desc(mapq)) %>%
  mutate(count=cumsum(count)) %>%
  ungroup()

input_csv %>%
  group_by(simulator, mapper) %>%
  mutate(total=sum(count)) %>%
  filter(unmapped==0) %>%
  select(-unmapped) %>%
  group_by(simulator, mapper, mapq, total) %>%
  summarize(count=sum(count)) %>%
  group_by(simulator, mapper) %>%
  arrange(desc(mapq)) %>%
  mutate(count=cumsum(count)) %>%
  ungroup() %>%
  rbind(input_real) %>%
  rename(Simulator=simulator,Mapper=mapper) %>%
  arrange(desc(Mapper)) %>% # Plot bowtie2 over bwa
  ggplot(aes(x=mapq, y=count*100/total, color=Mapper, shape=Simulator)) +
    geom_line(size=4) +
    geom_point(size=8) +
    scale_shape_manual(values=c(17, 19, 15)) +
    scale_color_manual(values=c("#009E73","#0072B2")) +
    scale_x_reverse() +
    xlab("Mapping quality") +
    ylab("% mapped") +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(0.99, 0.01),
          legend.justification = c(1, 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
          {if(nolegend)theme(legend.position = "None")}

ggsave(paste0(args[3],"_mapq.pdf"), width=297, height=210, units="mm")