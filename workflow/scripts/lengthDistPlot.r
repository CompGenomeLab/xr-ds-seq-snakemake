#### Libraries ####

library("ggplot2")
library(argparser)

######## Arguments ##########
p <- arg_parser("producing length distribution plot")
p <- add_argument(p, "-i", help="input")
p <- add_argument(p, "-o", help="output")
p <- add_argument(p, "-s", help="sample name written in the caption")
# Parse the command line arguments
argv <- parse_args(p)

#### Rearrange ####

# scientific notation
options(scipen=999)
  
d <- read.delim(argv$i, header = FALSE)

colnames(d) <- c("oligomer_length", "counts")

#### Plot ####

p <- ggplot(d, aes(x = oligomer_length, y = counts, 
                              fill = counts)) + 
  geom_bar(stat = "identity") +
  xlim(15, 35) +
  xlab("Oligomer Length") + ylab("Counts") +
  labs(title="Length Distribution of the Reads", 
       subtitle="Bar Chart",
       caption = argv$s)

# theme

p <- p + theme_light() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, vjust = 0.6),
        axis.text.y = element_text(size = 14, vjust = 0.1),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14, angle = 360),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        plot.caption = element_text(size = 14, 
                                    face = "italic"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 16, face = "bold"))

# save
ggsave(argv$o, width = 10, height = 8)

