library("alakazam")
library("ggplot2")
library("dplyr")
library('tidyverse')

#read metadata
meta <- read.csv("path/to/file/with/metadata.csv")
meta$ID <- gsub(" ","", meta$ID)

meta <- meta[((meta$CITY == "BH") | ((meta$CITY == "GV") & (meta$CLINICAL.CLASSIFICATION == "Mild"))),]
#list report file from clonotyping
clonotyped_report <- list.files(path="path/to/folder/with/tsv_YClon_clonotyped_files/", 
                                recursive = TRUE,
                                all.files = TRUE, 
                                pattern="clonotyped_report.tsv")

clonotyped <- list.files(path="path/to/folder/with/tsv_YClon_clonotyped_files/", 
                         recursive = TRUE,
                         all.files = TRUE, 
                         pattern="clonotyped.tsv")


unique_seq <- 0
names <- c("ID","Unique_sequences","Clones")
overall_char <- data.frame(matrix(ncol = length(names), nrow = 0))
colnames(overall_char) <- names
ct<-0
cdr3_len <- data.frame(ID="",junction_length="")
for(i in 1:length(clonotyped)){
  if(!strsplit(clonotyped[i],"/")[[1]][2] %in% unique(meta$ID)){
    next
  }
  print(strsplit(clonotyped[i],"/")[[1]][2])
  df <- read.csv(paste("path/to/folder/with/tsv_YClon_clonotyped_files/",clonotyped[i],sep=""),
                 sep="\t",)
  df <- df[df$productive==TRUE,]
  df <- df[!duplicated(df$sequence),]
  tmp <- data.frame(ID=strsplit(clonotyped[i],"/")[[1]][2],junction_length=df$junction_length)
  cdr3_len <- rbind(cdr3_len,tmp)

meta_len <- merge(cdr3_len,meta, on="ID")

meta_len$junction_length <- as.numeric(meta_len$junction_length)
meta_len$junction_length <- (meta_len$junction_length/3)-2

ggplot(meta_len[meta_len$CITY=="BH",]) +
  geom_density(mapping = aes(x=junction_length, color=CLINICAL.CLASSIFICATION),
               position = "identity",
               adjust = 4,
               linewidth=2) +
  theme_bw() +
  xlab("CDR3 Length AA") +
  ylab("Frequency") +
  labs(colour="Disease Outcome")
  

ggplot(meta_len[meta_len$CLINICAL.CLASSIFICATION=="Mild",]) +
  geom_density(mapping = aes(x=junction_length, color=CITY),
               position = "identity",
               adjust = 4,
               linewidth=2) +
  theme_bw() +
  xlab("CDR3 Length AA") +
  ylab("Frequency") +
  labs(colour="CITY") +
  scale_color_manual(values=c("#9999CC", "#66CC99")) +
  ylim(0,0.13)



