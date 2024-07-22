library("alakazam")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("iNEXT")
library('tidyverse')
library('vegan')
library('rstatix')
library('ggprism')
library("dbplyr")
library("tidyverse")


df <- read.table("path/to/tsv_YClon_clonotyped_file.tsv",
                    sep = "\t",
                    head = TRUE)

covid <- covid %>%
  select(clone_id)

covid <- covid %>% 
  add_column( clone_seq_count=NA )


df <- df %>%
  select(clone_id)

df <- df %>% 
  add_column( clone_seq_count=NA )


df <- covid[!duplicated(covid[ , c("clone_id")]),]
print("Calculating rarefaction... This will take some minutes")
x <- iNEXT(df$clone_seq_count)


p <-ggiNEXT(x,color.var="Order.q") +
  xlab("Number of sequences") +
  ylab("Clone diversity") +
  annotate("text", x = 3000, y =-3000, label = paste("Coverage: ",coverage), size=8) +
  coord_cartesian(ylim=c(-0,30000),clip="off")


print("Saving the plot!")
png(str_replace(args[1],".tsv","_RAREFACTION.png"),width = 1200, height = 1200,)
print(p)
dev.off()