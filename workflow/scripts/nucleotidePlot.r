#### Libraries ####

library(ggplot2)
library(argparser)
library("tidyr")
library(dplyr)
library(futile.logger)

######## Arguments ##########
p <- arg_parser("producing di/nucleotide contents")
p <- add_argument(p, "-i", help="input")
p <- add_argument(p, "-o", help="output")
p <- add_argument(p, "-l", help="log file")
p <- add_argument(p, "-k", help="kmer value of input, 1 (nucleotide) or 2 (dinucleotide)")
p <- add_argument(p, "-s", help="sample name written in the caption")
p <- add_argument(p, "-f", help="filter (only for dinucleotides). Ex: 'CC','CT','TC','TT'")
# Parse the command line arguments
argv <- parse_args(p)

# log file
flog.appender(appender.file(argv$l))

#### Rearrange ####

flog.info("Reading file...")  
nuc_table <- read.table(argv$i, header = TRUE)

flog.info("Renaming columns...")  
x_order <- c()

if (argv$k == "1"){
  for (i in 2:ncol(nuc_table)) {
    
    colnames(nuc_table)[i] <- c(paste(i-1))
    
    x_order <- c(x_order, paste(i-1)) 
    
  }
  
  #title for the plot
  mytitle <- "Nucleotide Content of Oligomers"

} else if (argv$k == "2"){
  for (i in 2:ncol(nuc_table)) {
    
    colnames(nuc_table)[i] <- c(paste(i - 1, "-", i, sep = ""))
    
    x_order <- c(x_order, paste(i - 1, "-", i, sep = "")) 
    
  }

  #title for the plot
  mytitle <- "Dinucleotide Content of Oligomers"
}

# reorganize
dt_organized <- nuc_table %>% gather(Position, count, 
                                              2:ncol(nuc_table))

colnames(dt_organized) <- c("nucleotides", "positions", "counts")

flog.info("Calculating the percentages...")  
dt_organized$freq <- 100*dt_organized$counts/sum(nuc_table[2])

if (argv$k == "2"){
  
  new_dt <- data.frame()
  filt <- unlist(strsplit(argv$f, ","))
  flog.info("Motifs: %s", list(filt))  
  for (i in filt){
    new_dt <- rbind(new_dt, filter(dt_organized, nucleotides == i))
    
  }
  dt_organized <- unique(new_dt)
}

flog.info("Ordering columns...")  
dt_organized$positions = factor(dt_organized$positions, 
                                levels = x_order)

#### Plot ####

flog.info("Plotting...")  
p <- ggplot(dt_organized, aes(x = positions, y = freq, 
                              fill = nucleotides)) + 
  geom_bar(stat = "identity") +
  ylim(0,101) +
  ylab("Frequency (%)") + xlab("Position in oligomers") +
  labs(title = mytitle, 
       subtitle = "Categorywise Bar Chart", 
       caption = argv$s)

# theme 
p <- p + theme_light() +
  theme(axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 14, vjust = 0.6, angle = 65),
        axis.text.y = element_text(size = 14, vjust = 0.1),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14, angle = 360),
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 14),
        plot.caption = element_text(size = 14, 
                                    face = "italic"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 16, face = "bold"))

ggsave(argv$o, width = 10, height = 8)

flog.info("Saved.")  