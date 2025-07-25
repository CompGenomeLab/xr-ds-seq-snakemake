#### Libraries ####

library("ggplot2")
library(argparser)
library(logger)
set.seed(1)

######## Arguments ##########
p <- arg_parser("producing the length distribution plot")
p <- add_argument(p, "-i", help = "input")
p <- add_argument(p, "-o", help = "output")
p <- add_argument(p, "-l", help = "log file")
p <- add_argument(
  p, "-s",
  help = "sample name written in the caption"
)
# Parse the command line arguments
argv <- parse_args(p)

# log file
log_appender(appender_file(argv$l))

#### Rearrange ####

# scientific notation
options(scipen = 999)

log_info("Reading file...")
d <- read.delim(argv$i, header = FALSE)

colnames(d) <- c("oligomer_length", "counts")

#### Plot ####

most_occ_oligo <- d[which.max(d[, 2]), 1]

log_info("Plotting...")
p <- ggplot(d, aes(
  x = oligomer_length, y = counts, fill = counts
)) +
  geom_bar(stat = "identity") +
  xlim(most_occ_oligo - 10, most_occ_oligo + 10) +
  xlab("Oligomer Length") +
  ylab("Counts") +
  labs(
    title = "Length Distribution of the Reads",
    subtitle = "Bar Chart",
    caption = argv$s
  )

# theme
p <- p + theme_light() +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 14, vjust = 0.6),
    axis.text.y = element_text(size = 14, vjust = 0.1),
    strip.text.x = element_text(size = 14),
    strip.text.y = element_text(size = 14, angle = 360),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    plot.caption = element_text(size = 14, face = "italic"),
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 16, face = "bold")
  )

ggsave(argv$o, width = 10, height = 8)

log_info("Saved.")
