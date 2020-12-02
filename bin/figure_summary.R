suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(purrr))
suppressMessages(library(stringr))
suppressMessages(library(patchwork))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))

args = commandArgs(trailingOnly=TRUE)

if(length(args)<3){
  stop("Csv file or one of the two output files missing", call.=FALSE)
} else if(length(args)>3){
  stop("Too many parameters. Only csv file and two output files required.", call.=FALSE)
}

# Input is strand of read, type is strand of fragment
summary_csv <- read_csv(args[1], col_types = cols())

order <- summary_csv %>%
  select(-Type) %>%
  gather(Simulator, level, -Comparison, -Dataset) %>%
  group_by(Simulator) %>%
  summarize(level=mean(level, na.rm = TRUE)) %>%
  arrange(level) %>%
  select(Simulator) %>%
  pull

text_size = 24
tick_text_size <- 24
data <- summary_csv %>%
  mutate(Comparison=factor(Comparison, levels=summary_csv %>% select(Comparison) %>% unique() %>% pull())) %>%
  mutate(Dataset=factor(Dataset, levels=summary_csv %>% select(Dataset) %>% unique() %>% pull())) %>%
  gather(Simulator, level, -Comparison, -Dataset, -Type) %>%
  mutate(Simulator=factor(Simulator, levels=order)) %>%
  mutate(level = recode(as.factor(level), "1"="very good", "2"="good", "3"="intermediate", "4"="poor", "5"="very poor", "NA"="NA")) %>%
  mutate(x = (as.integer(Dataset)-1)%%3, y = 3-((as.integer(Dataset)-1)%/%3)) %>%
  mutate(category=if_else(Comparison %in% c("K-mer spectrum", "Assembly N50", "Mapping qualities"), 2, if_else(str_detect(Comparison, "Speed|Memory"), 3, 1)))

get_plot <- function(cdata, cat){
  cdata %>%
  filter(category == cat) %>%
  ggplot(aes(Simulator, Comparison, fill = level)) +
  geom_tile(aes(x, y, fill = level)) +
  geom_text(aes(x, y, label=tile_text), size=9) +
  scale_fill_manual(values=c("#488BC2","#7FB972","#D9AD3C","#E6642C","#D92120"), na.value="#BBBBBB", drop=FALSE) +
  scale_x_discrete(expand = c(0.00,0.00)) +
  scale_y_discrete(expand = c(0.00,0.00)) +
  facet_grid( Comparison ~ Simulator, switch="y" ) +
  theme_bw() +
  theme(  axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text.x = element_text( size = text_size),
          strip.text.y = element_text( size = text_size, angle = 180, hjust=1),
          strip.background = element_blank(),
          panel.spacing=unit(0.1, "lines"),
          legend.title=element_blank(), 
          legend.text=element_text(size=text_size)) +
          {if(cat != 1) theme(strip.text.x = element_blank())}
}

plots <- lapply(1:3, function(cat){
  data %>%
    filter(Type == "Dataset replication") %>%
    mutate(tile_text = if_else(as.integer(Comparison)==1 & as.integer(Simulator)==max(as.integer(Simulator)), LETTERS[as.integer(Dataset)], "")) %>% # tolower(LETTERS[as.integer(Dataset)])
    get_plot(cat)
  })

datasets <- summary_csv %>%
  filter(Type == "Dataset replication") %>%
  select(Dataset) %>%
  unique() %>%
  pull()

plots[[1]] + plots[[2]] + plots[[3]] +
#  inset_element(
#    tableGrob(data.frame(c("a)","b)","c)","d)","e)","f)","g)","h)","i)","j)","k)","l)"), datasets), 
#              theme=ttheme_minimal(core = list(fg_params = list(hjust=0, x=0, fontsize=tick_text_size))),
#              rows=NULL, cols=NULL),
#    left=1.255, bottom=1.9, right=1.555, top=3.4) +
  plot_layout(heights = c(7/14, 3/14, 4/14), guides = "collect") &
  theme( plot.margin = unit(c(0.4,0,0,0), "lines"),
         legend.position = "none")
#         legend.justification = "bottom",
#         legend.margin = margin(t=0, r=120, b=0, l=0),
#         legend.key.height = unit(2, "lines"))

ggsave(args[2], width=297, height=420, units="mm")

plots <- lapply(1:2, function(cat){
  data %>%
    filter(Type == "Cross species simulation") %>%
    mutate(tile_text = if_else(as.integer(Comparison)==1 & as.integer(Simulator)==max(as.integer(Simulator)), LETTERS[as.integer(Dataset)], "")) %>%
    get_plot(cat)
})

table <- summary_csv %>%
  filter(Type == "Cross species simulation") %>%
  select(Dataset) %>%
  unique() %>%
  mutate(let = c("a)","b)","c)","d)","e)","f)")) %>%
  separate(Dataset, c("dataset", "profile"), sep=" <- ") %>%
  mutate(profile = paste0("<- ", profile)) %>%
  gather(key, set, -let) %>%
  arrange(let) %>%
  mutate(let=if_else(key=="profile", "", let)) %>%
  select(-key)

plots[[1]] + plots[[2]] +
  #inset_element(
  #  tableGrob( table, 
  #            theme=ttheme_minimal(core = list(fg_params = list(hjust=0, x=0, fontsize=tick_text_size))),
  #            rows=NULL, cols=NULL),
  #  left=1.18, bottom=1.1, right=1.48, top=2.6) +
  plot_layout(heights = c(48/80, 32/80), guides = "collect") &
  theme( plot.margin = unit(c(0.4,0,0,0), "lines"),
         #legend.justification = "bottom",
         #legend.margin = margin(t=0, r=120, b=0, l=0),
         legend.key.height = unit(2, "lines"))

ggsave(args[3], width=297, height=420, units="mm")







