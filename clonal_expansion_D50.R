library(Biostrings) 
library(kmer)
library(igraph)
library("dplyr")
library("tibble")
library('tidyverse')
library(ape)
library(multivariance)
library(plyr)
library(data.table)
library("qgraph")

library("colorspace")

max_min <- function(x){
  if(x*0.2<1){
    return(0.3)
  }else if(x*0.2>10){
    return(5)
  }else{
    return(1)
  }
}
meta <- read.csv("path/to/file/with/metadata.csv")
meta$ID <- gsub(" ","", meta$ID)
meta <- meta[((meta$CITY=="BH") & (meta$CLINICAL.CLASSIFICATION=="Mild")),]

df <- read.csv(paste("path/to/folder/with/tsv_YPub_clonotyped_files/","YPub_input_YClon_clonotyped.tsv",sep="/"),
               sep="\t") 
df <- df[!(duplicated(df[c("sequence","sample_ID")])),]
df <- df %>% drop_na(clone_id)
a <- df %>%
  dplyr::group_by(clone_id) %>%
  dplyr::summarize(publicity = n_distinct(sample_ID))

df <- left_join(df,a,by="clone_id")

meta$sample_ID <- gsub(" ","", meta$ID)
df <- df[df$sample_ID %in% meta$sample_ID,]

for(i in unique(df$sample_ID)){
  if(i!="ID117"){
    next
  }
  sample <- i
  print(i)
  df_red <- df[df$sample_ID==sample,]
  
  df_red <- df_red %>% dplyr::group_by(clone_id) %>% 
    dplyr::mutate(clone_seq_count = n()) %>%
    dplyr::ungroup()
  
  unique_sum <- df_red %>%
    dplyr::distinct(clone_id, clone_seq_count, .keep_all = TRUE) %>%
    dplyr::group_by(clone_id) %>%
    dplyr::summarise(Sum_Unique_Values = sum(clone_seq_count)) %>%
    dplyr::ungroup() %>%
    arrange(desc(Sum_Unique_Values)) %>%
    dplyr::mutate(cum_sum = cumsum(Sum_Unique_Values))
     
  full_rep = max(unique_sum$cum_sum)
  d50_clones <- which.min(abs(unique_sum$cum_sum - (full_rep/2)))
  d50_clones <- unique_sum[1:d50_clones, ]
  df_red <- df_red[df_red$clone_id %in% d50_clones$clone_id,]
  df_red <- df_red[df_red$clone_seq_count>5,]
  
  cols <-c("#EFCFE5")
  full_size <- list()
  pub_color <- list()
  final_mst <- 'A' 
  
  for(x in unique(df_red$clone_id)){
    clone_dist_mat <- df_red[df_red$clone_id==x,]
    clone_dist_mat <- clone_dist_mat %>% drop_na(sequence)
    if(nrow(clone_dist_mat) <=1){
      next
    }
    dna_string <- DNAStringSet(clone_dist_mat$sequence)
    names(dna_string)<- clone_dist_mat$sequence_id
    sequences <- as.DNAbin(dna_string)
    ktest <- kdistance(sequences, k = 4, gap=".")
    dist <- as.matrix(ktest)
    unique_nodes <- !duplicated(dist, fromLast = TRUE)
    clone_dist_mat$color <- ifelse(clone_dist_mat$publicity > 1, darken(cols, clone_dist_mat$publicity*0.1), "blue")
    pub_color <- append(pub_color,clone_dist_mat$color)
    size <- rowSums(data.frame(dist)  == 0)
    size <- unlist(lapply(size, max_min))
    full_size <- append(full_size,size)
    
    # dist_matrix_tmp <- dist_matrix
    dist_matrix <- dist[unique_nodes, unique_nodes]
    
    tree_graph <- graph.adjacency(dist_matrix, mode = "undirected", weighted = TRUE)
    mst <- minimum.spanning.tree(tree_graph)
    if(typeof(final_mst)=="character"){
      attrs <- get.data.frame(mst, what= c("vertices")) %>% unique()
      el <- get.data.frame(mst, what= c( "edges"))
      final_mst <- mst
    }else{
      attrs <- rbind(attrs,get.data.frame(mst, what= c("vertices"))) %>% unique()
      el <- rbind(el,get.data.frame(mst, what= c( "edges")))
      final_mst <- mst
    }
    
    
  }
  g3 <- graph_from_data_frame(el, directed = FALSE, vertices = attrs)
  e <- get.edgelist(g3, names=FALSE)
  l <- qgraph.layout.fruchtermanreingold(e,
                                         vcount=vcount(g3),
                                         repulse.rad=vcount(g3)**4,
                                         niter=1000)  
  
  png(filename = paste('path/to/folder/with/samples_results/expansion_D50_images/',i,'clonal_exp.png',sep=""),width = 6500, height = 5767,
      units = "px")
  plot(g3,layout = l,
       vertex.size = unlist(full_size),
       vertex.color = unlist(pub_color),
       vertex.label = NA,
       edge.label=NA)
  dev.off() 



}