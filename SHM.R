library("ggplot2")
library("dplyr")
library('tidyverse')


#read metadata
meta <- read.csv("path/to/file/with/metadata.csv")
meta$ID <- gsub(" ","", meta$ID)


  
  ########### SHM ########
  library(shazam)
  library(alakazam)
  library(dplyr)
  library('rstatix')
  library(ggpubr)
  library('ggprism')
  
  clonotyped <- list.files(path="path/to/folder/with/tsv_YClon_clonotyped_files/", 
                           recursive = TRUE,
                           all.files = TRUE, 
                           pattern="db-pass_YClon_clonotyped.tsv")
  
  
  SHM_all <- "A"

  
  for (x in clonotyped){
    if (unlist(strsplit(x,"/"))[2] %in% unique(SHM_all$ID)){
      next
    }else{
      print(unlist(strsplit(x,"/"))[2])
      df <- read.csv(paste("path/to/folder/with/tsv_YClon_clonotyped_files",x,sep="/"),
                     sep="\t")
      
      df <- df[df$clone_seq_count>10,]
      
      df <- df %>%
        select(sequence_id,clone_id,sampleID,sequence_alignment,germline_alignment)
      
      db_obs <- observedMutations(df, sequenceColumn="sequence_alignment",
                                  germlineColumn="germline_alignment",
                                  regionDefinition=NULL,
                                  frequency=TRUE,
                                  cloneColumn = "clone_id",
                                  nproc=4)
      db_obs["ID"] <- unlist(strsplit(x,"/"))[2]
      
      if(typeof(SHM_all)=="character"){
        SHM_all <- db_obs
      }else{
        SHM_all <- rbind(db_obs,SHM_all)
      }
    }
    
    
    
    
  }

final_result <- left_join(SHM_all,meta, by="ID")

clonotyped_report <- list.files(path="path/to/folder/with/tsv_YClon_clonotyped_files/", 
                                recursive = TRUE,
                                all.files = TRUE, 
                                pattern="_clonotyped_report.tsv")


mean_mu_r <- final_result %>%
  group_by(ID,clone_id) %>%
  summarize(mu_r_mean = mean(mu_freq_seq_r), .groups = 'drop')


mean_mu_r <- left_join(mean_mu_r,meta, by="ID")
mean_p_val <- rstatix::wilcox_test(mean_mu_r[mean_mu_r$CITY=="BH",], mu_r_mean ~ CLINICAL.CLASSIFICATION) %>%
  rstatix::add_xy_position()
ggplot(mean_mu_r[mean_mu_r$CITY=="BH",], aes(x=CLINICAL.CLASSIFICATION,y=mu_r_mean)) +
  geom_boxplot(aes(fill=CLINICAL.CLASSIFICATION)) +
  ylab("Nonsynonymous Mutational Frequency") +
  xlab("Disease Outcome") +
  labs(fill="Disease Outcome", title = "BH only") +
  add_pvalue(mean_p_val[mean_p_val$p<0.05,], label = "p.adj.signif", step.increase = 0.05) +
  theme_bw()



mean_mu_s <- final_result %>%
  group_by(ID,clone_id) %>%
  summarize(mu_s_mean = mean(mu_freq_seq_s), .groups = 'drop')


mean_mu_s <- left_join(mean_mu_s,meta, by="ID")
mean_p_val <- rstatix::wilcox_test(mean_mu_s[mean_mu_s$CITY=="BH",], mu_s_mean ~ CLINICAL.CLASSIFICATION) %>%
  rstatix::add_xy_position()
ggplot(mean_mu_s[mean_mu_s$CITY=="BH",], aes(x=CLINICAL.CLASSIFICATION,y=mu_s_mean)) +
  geom_boxplot(aes(fill=CLINICAL.CLASSIFICATION)) +
  ylab("Synonymous Mutational Frequency") +
  xlab("Disease Outcome") +
  labs(fill="Disease Outcome", title = "BH only") +
  add_pvalue(mean_p_val[mean_p_val$p<0.05,], label = "p.adj.signif", step.increase = 0.05) +
  theme_bw()


mean_mu <- final_result %>%
  group_by(ID,clone_id) %>%
  summarize(mu_mean = mean(mu_freq_seq_s+mu_freq_seq_r), .groups = 'drop')


mean_mu <- left_join(mean_mu,meta, by="ID")
mean_p_val <- rstatix::wilcox_test(mean_mu[mean_mu$CITY=="BH",], mu_mean ~ CLINICAL.CLASSIFICATION) %>%
  rstatix::add_xy_position()
ggplot(mean_mu[mean_mu$CITY=="BH",], aes(x=CLINICAL.CLASSIFICATION,y=mu_mean)) +
  geom_boxplot(aes(fill=CLINICAL.CLASSIFICATION)) +
  ylab("Overall Mutational Frequency") +
  xlab("Disease Outcome") +
  labs(fill="Disease Outcome", title = "BH only") +
  add_pvalue(mean_p_val[mean_p_val$p<0.05,], label = "p.adj.signif", step.increase = 0.05) +
  theme_bw()
