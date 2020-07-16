suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library("viridis"))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<4){
  stop("Load directory and three output files must be supplied", call.=FALSE)
} else if(length(args)>4){
  stop("Only load directory and three output file must be supplied", call.=FALSE)
}

input_path = args[1]

ipf_comp <- read_csv(paste0(input_path,"SRR490124_seg0_a_er0_base_quality_preceding_quality_sequence_quality_position.csv.gz"), col_types = cols()) %>%
  full_join( read_csv(paste0(input_path,"SRR490124_estimate_seg0_a_er0_base_quality_preceding_quality_sequence_quality_position.csv.gz"), col_types = cols()), by = c("base_quality", "previous_quality", "sequence_quality", "position") ) %>%
  mutate( count = replace_na(count, 0) ) %>%
  group_by(previous_quality, sequence_quality, position) %>%
  mutate(total_counts = sum(count)) %>%
  ungroup %>%
  filter(0 < total_counts) %>%
  mutate( likelihood = replace_na(likelihood, 0.0) )

# Check interesting variable selections
#count_overview <- ipf_comp %>% 
#  group_by(previous_quality, sequence_quality, position) %>%
#  summarize(total_counts=first(total_counts)) %>%
#  ungroup() %>%
#  arrange(desc(total_counts))
  
#count_overview %>% head(100) %>% View()
#count_overview %>%
#  filter(10000 <= total_counts) %>%
#  tail(1)
#count_overview %>%
#  filter(10 <= total_counts) %>%
#  tail(100) %>% View()

ipf_comp_plot <- ipf_comp %>%
  group_by(previous_quality, sequence_quality, position) %>%
  mutate(total_likelihood = sum(likelihood)) %>%
  ungroup %>%
  mutate(Estimation=likelihood/total_likelihood*total_counts) %>%
  rename(Observation=count) %>%
  gather(Type,Counts,Observation,Estimation)

text_size <- 20
tick_text_size <- 16

ipf_comp_plot %>%
  filter(   (1 == previous_quality & 37 == sequence_quality & 0 == position) |
            (39 == previous_quality & 37 == sequence_quality & 50 == position) |
            (34 == previous_quality & 34 == sequence_quality & 90 == position) |
            (29 == previous_quality & 24 == sequence_quality & 99 == position) ) %>%
  mutate(group = paste0("N=",total_counts,ifelse(0==position,"",paste0(", PQ=",previous_quality)),", SQ=",sequence_quality,", Pos=",position)) %>%
  ggplot(aes(x=base_quality, y=Counts, colour=Type)) +
    geom_point(size=3) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), limits=c(1e-2, NA)) +
    scale_color_manual(values=c("#488BC2","#E6642C")) +
    xlab("Base Quality") +
    facet_wrap(. ~ group) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(0.01, 0.99),
          legend.justification = c(0, 1))

ggsave(args[2], width=297, height=320, units="mm")

ipf_comp_plot2 <- ipf_comp %>%
  group_by(previous_quality, sequence_quality, position) %>%
  mutate(total_likelihood = sum(likelihood)) %>%
  ungroup %>%
  mutate(prob=count/total_counts, est_prob=likelihood/total_likelihood) %>%
  mutate(diff=est_prob-prob, rel_diff=diff/est_prob)

text_size <- 24
tick_text_size <- 20

ipf_comp_plot2 %>%
  filter(0 < count) %>%
  mutate( Group = ifelse(count < 10, "[1,10[", ifelse(count < 100, "[10,10^2[", ifelse(count < 1000, "[10^2,10^3[", ifelse(count < 10000, "[10^3,10^4[", ifelse(count < 100000, "[10^4,10^5[", ">= 10^5"))))) ) %>%
  group_by(Group) %>% mutate(num_rows = n(), group_counts=sum(count)) %>% ungroup() %>%
  mutate( Group = sprintf("%s\nN=%d\n%.2f%%", Group, num_rows, group_counts/sum(count)*100) ) %>%
  mutate( type = ifelse( 0 <= rel_diff, "Overestimated", "Underestimated") ) %>%
  ggplot( aes(x=Group, y=abs(rel_diff), fill=Group) ) +
    geom_boxplot(na.rm = TRUE) +
    scale_y_log10(limits=c(10e-8,100), breaks = trans_breaks("log10", function(x) 10^x)) +
    scale_fill_viridis(discrete = TRUE) +
    labs(y=expression(frac(paste("|",Estimated-Observed,"|"),Estimated)),x="Group") +
    facet_grid(. ~ type) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.key.size = unit(1, "cm"),
          legend.key.height=unit(3,"cm"))

ggsave(args[3], width=297, height=210, units="mm")

ipf_comp_plot2 %>%
  arrange(abs(diff)) %>%
  mutate(cumcount=round(cumsum(count)/sum(count)*200)/2) %>%
  mutate( Group = ifelse(count < 10, "[0,10[", ifelse(count < 100, "[10,10^2[", ifelse(count < 1000, "[10^2,10^3[", ifelse(count < 10000, "[10^3,10^4[", ifelse(count < 100000, "[10^4,10^5[", ">= 10^5"))))) ) %>%
  mutate( type = ifelse( 0 <= diff, "Overestimated", "Underestimated") ) %>%
  group_by(type,Group,cumcount) %>%
  summarize(diff=max(abs(diff))) %>%
  ungroup() %>%
  ggplot( aes(x=diff, y=cumcount, colour=Group) ) +
    geom_point(size=6) +
    scale_x_continuous(breaks=1:10/10) +
    scale_color_viridis(discrete = TRUE) +
    labs(x=expression("|estimated probability - observed fraction|") , y="% observed counts") +
    facet_grid(type ~ .) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=text_size),
          legend.position = c(0.99, 0.01),
          legend.justification = c(1, 0))

ggsave(args[4], width=297, height=210, units="mm")
