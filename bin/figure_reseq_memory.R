suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<2){
  stop("Input and output file must be supplied", call.=FALSE)
} else if(length(args)>2){
  stop("Only input and output file must be supplied", call.=FALSE)
}

text_size <- 22
tick_text_size <- 18

mem_csv <- read_csv(args[1], col_types = cols())

mem_csv %>%
  filter(Simulator == "ReSeq") %>%
  mutate(max_memory = max_memory/1024) %>%
  ggplot(aes(x=threads, y=max_memory, color=Simulator)) +
    geom_line(na.rm=TRUE, size=4) +
    geom_point(na.rm=TRUE, size=8) +
    scale_color_manual(values=c("#488BC2")) +
    scale_x_continuous(breaks=c(1,4,8,16,32,48,64)) +
    scale_y_continuous(limits=c(0,256), breaks=c(0,16,32,64,128,192,256)) +
    xlab("Threads") +
    ylab("Max. memory [GB]") +
    ggtitle("Training\nHs-HiX-TruSeq (3.3Gb, 40x)") +
    theme_bw() +
    theme(axis.text.x = element_text( size = tick_text_size, angle = 90, vjust=0.5),
          axis.text.y = element_text( size = tick_text_size),
          axis.title.x = element_text( size = text_size),
          axis.title.y = element_text( size = text_size),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size),
          plot.title = element_text( size= text_size, hjust=0.5),
          legend.title=element_text(size=text_size), 
          legend.text=element_text(size=tick_text_size),
          legend.justification = c(0, 1),
          legend.position = c(0.01, 0.99),
          panel.grid.minor = element_blank())

ggsave(args[2], width=297, height=210, units="mm")