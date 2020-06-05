suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<3){
  stop("Input file, dataset and output file must be supplied", call.=FALSE)
} else if(length(args)>3){
  stop("Only input file, dataset and output file must be supplied", call.=FALSE)
}

map_stats <- read_csv(args[1], col_types = cols())
nolegend <- grepl("nolegend", args[3], fixed=TRUE)

text_size <- 28
tick_text_size <- 24

sample <- recode(args[2], "SRR490124"="Ec-Hi2000-TruSeq", "SRR3191692"="Ec-Hi2500-TruSeq", "S5L001"="Ec-Hi4000-Nextera")
map_stats %>%
  gather(key,Value,-Variable,-Simulator) %>%
  separate(key, c("Dataset","Mapper")) %>%
  mutate(group=ifelse(Variable %in% c("Unmapped pairs","Single reads","Mapped pairs"), "Read", "Cigar")) %>%
  group_by(Dataset, Simulator, Mapper, group) %>%
  mutate(total = sum(Value), percent=Value/total*100) %>%
  ungroup() %>%
  filter(!(Variable %in% c("Mapped pairs","Matches")), Dataset == args[2], percent != 0.0) %>%
  ggplot(aes(x=Variable, y=percent, color=Mapper, shape=Simulator)) +
    geom_point(na.rm=TRUE, size=8, position=position_dodge(width = 0.6)) +
    scale_y_log10() +
    scale_shape_manual(values=c(17, 19, 15)) +
    scale_color_manual(values=c("#009E73","#0072B2")) +
    facet_grid(. ~ group, scales = "free") +
    xlab("") +
    ylab("Percent of total") +
    ggtitle(sample) +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size, angle = 20, vjust=0.5),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          plot.title = element_text( size= text_size, hjust=0.5),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.position = c(0.01, 0.99),
          legend.justification = c(0, 1),
          panel.grid.minor = element_blank()) +
    {if(nolegend)theme(legend.position = "None")}


ggsave(args[3], width=297, height=210, units="mm")