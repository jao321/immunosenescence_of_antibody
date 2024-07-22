library("alakazam")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("iNEXT")
library('tidyverse')
library('vegan')
library('rstatix')
library('ggprism')
library("tibble")
library("dplyr")
# library("poppr")
# library("vignette")


#DIFERENÇA ESTATÍSTICA ENTRE LEVE NAS 3 CIDADES*******

clonotyped <- list.files(path="path/to/folder/with/tsv_YClon_clonotyped_files/", 
                recursive = TRUE,
                all.files = TRUE, 
                pattern="clonotyped_report.tsv")

diversities <- data.frame()

diversities <- diversities %>% 
  add_column( ID=NA ) %>%
  add_column( shannon=NA ) %>%
  add_column( simpson=NA )


meta <- read.csv("path/to/file/with/metadata.csv")
meta$ID <- gsub(" ","", meta$ID)

for(i in clonotyped){
  df <- read.csv(paste("path/to/folder/with/tsv_YClon_clonotyped_files/",i,sep=""),
                 sep="\t",
                 header = FALSE)
  colnames(df)=c("sequence_id","seq_count","most_common_cdr3","clone_id")
  df <- df[order(df$seq_count, decreasing = TRUE),]
  row.names(df) <- NULL
  df <- df %>%
    add_column(cumsum=cumsum(df$seq_count))
  full_rep = sum(df$seq_count)
  d50_clones <- as.numeric(which.min(abs((full_rep/2)-df$cumsum)))
  cat(strsplit(i,"/")[[1]][2],",",d50_clones,"\n")
  df50 <- df[1:d50_clones,]
  tmp <- data.frame(strsplit(i,"/")[[1]][2],diversity(df50$seq_count, "shannon"),simpson.unb(df50$seq_count))
  colnames(tmp) <- c("ID","shannon","simpson")
  diversities <- rbind(diversities, tmp)
}


diversities <- left_join(diversities,meta%>%select(ID,AGE,CLINICAL.CLASSIFICATION,AGE_GROUP,CITY), by = "ID")

diversities[diversities$CLINICAL.CLASSIFICATION=="Control",]$CITY="Control"
diversities[diversities$CLINICAL.CLASSIFICATION=="Control",]$AGE_GROUP="Control"

ggplot(diversities, aes(color=CLINICAL.CLASSIFICATION,x=AGE,y=shannon,)) +
  geom_count(size=5) +
  scale_color_manual(values = c("#1B79A5", "#FD7701","#FF4E00","#8ea604")) +
  geom_smooth(method="lm",se=FALSE)

diversities_city <- diversities[(((diversities$CITY=="BH")|(diversities$CITY=="GV"))
                                & ((diversities$CLINICAL.CLASSIFICATION=="Mild"))),]
df_p_val <- rstatix::wilcox_test(diversities_city, shannon ~ CITY) %>%
  rstatix::add_xy_position()


annotation <- data.frame(
  x = c(1,2),
  y = c(1,1),
  label = c(paste("n=",
                  nrow(diversities[diversities$CITY=="BH",])),
            paste("n=",
                  nrow(diversities[diversities$CITY=="GV",])))
)

df_p_val <- diversities_city %>% rstatix::t_test(shannon ~ CITY) %>%
  rstatix::add_xy_position()

diversities$CITY <- factor(diversities$CITY, 
                           levels = c("Control","BH","GV","SP"))

result <- t.test(diversities_city[diversities_city$CITY=="GV",'shannon'],
       diversities_city[diversities_city$CITY=="BH",'shannon'])
diversities_city <- na.omit(diversities_city)
result <- diversities_city %>%
  t_test(shannon ~ CITY, var.equal = TRUE)


library(ggpubr)

ggplot(diversities_city, aes(x = CITY, y = shannon )) + 
  ylab("Shannon Diversity") +
  xlab("City") +
  labs(fill="City") +
  geom_boxplot(aes(fill=CITY)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  stat_pvalue_manual(
    result, 
    label = "p = {p}",
    tip.length = 0.05,
    hide.ns = TRUE
  ) +
  add_pvalue(df_p_val[df_p_val$p<0.06,], label = "p", step.increase = 0.03) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = alpha("white", 0)),  # Transparent background
    panel.grid.major = element_line(color = "black", size = 0.2),  # Black grid lines
    panel.grid.minor = element_blank(),  # No minor grid lines
    axis.line = element_line(color = "black"),
  )


annotation <- data.frame(
  x = c(1,2,3,4),
  y = c(0.7,0.7,0.7,0.7),
  label = c(paste("n=",
                  nrow(diversities[diversities$CITY=="Control",])),
            paste("n=",
                  nrow(diversities[diversities$CITY=="BH",])),
            paste("n=",
                  nrow(diversities[diversities$CITY=="GV",])),
            paste("n=",
                  nrow(diversities[diversities$CITY=="SP",])))
)
df_p_val <- rstatix::wilcox_test(diversities, simpson ~ CITY) %>%
  rstatix::add_xy_position()

ggplot(diversities, aes(x = CITY, y = simpson )) + 
  ylab("Simpson Diversity") +
  xlab("City") +
  labs(fill="City") +
  geom_boxplot(aes(fill=CITY)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" ) +
  add_pvalue(df_p_val[df_p_val$p<0.06,], label = "p", step.increase = 0.03)
  # add_pvalue(df_p_val[df_p_val$p.adj.signif!="ns",], 
  #            label = "p.adj.signif",
  #            y.position = 2,
  #            tip.length = -0.02) 
   



df_p_val <- rstatix::wilcox_test(diversities, shannon ~ AGE_GROUP) %>%
  rstatix::add_xy_position() 



diversities$AGE_GROUP <- factor(diversities$AGE_GROUP, 
                           levels = c("Control","Adult","Elderly"))



df_p_val <- rstatix::wilcox_test(diversities, shannon ~ AGE_GROUP) %>%
  rstatix::add_xy_position() 


annotation <- data.frame(
  x = c(1,2,3),
  y = c(1,1,1),
  label = c(paste("n=",
                  nrow(diversities[diversities$AGE_GROUP=="Control",])),
            paste("n=",
                  nrow(diversities[diversities$AGE_GROUP=="Adult",])),
            paste("n=",
                  nrow(diversities[diversities$AGE_GROUP=="Elderly",])))
)


ggplot(diversities, aes(x = AGE_GROUP, y = shannon )) + 
  ylab("Shannon Diversity") +
  xlab("Age Group") +
  labs(fill="Age Group") +
  geom_boxplot(aes(fill=AGE_GROUP)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" ) +
  add_pvalue(df_p_val[df_p_val$p<0.06,], label = "p", step.increase = 0.03)+
  add_pvalue(df_p_val[df_p_val$p.adj.signif!="ns",],
             label = "p.adj.signif",
             y.position = 2,
             tip.length = -0.02)




annotation <- data.frame(
  x = c(1,2,3),
  y = c(0.7,0.7,0.7),
  label = c(paste("n=",
                  nrow(diversities[diversities$AGE_GROUP=="Control",])),
            paste("n=",
                  nrow(diversities[diversities$AGE_GROUP=="Adult",])),
            paste("n=",
                  nrow(diversities[diversities$AGE_GROUP=="Elderly",])))
)


df_p_val <- rstatix::wilcox_test(diversities, simpson ~ AGE_GROUP) %>%
  rstatix::add_xy_position() 

ggplot(diversities, aes(x = AGE_GROUP, y = simpson )) + 
  ylab("Simpson Diversity") +
  xlab("Age Group") +
  labs(fill="Age Group") +
  geom_boxplot(aes(fill=AGE_GROUP)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" ) +
  add_pvalue(df_p_val[df_p_val$p<0.06,], label = "p",
             step.increase = 0.03)+
  add_pvalue(df_p_val[df_p_val$p.adj.signif!="ns",],
             label = "p.adj.signif",
             y.position = 0.75,
             tip.length = -0.02,
             step.increase = 0.03)

diversities[diversities$CLINICAL.CLASSIFICATION=="Moderate",]$CLINICAL.CLASSIFICATION="Hospitalized"
diversities[diversities$CLINICAL.CLASSIFICATION=="Severe",]$CLINICAL.CLASSIFICATION="Hospitalized"
diversities[is.na(diversities$CLINICAL.CLASSIFICATION),]$CLINICAL.CLASSIFICATION="Hospitalized"

df_p_val <- rstatix::wilcox_test(diversities, shannon ~ CLINICAL.CLASSIFICATION) %>%
  rstatix::add_xy_position()


diversities$CLINICAL.CLASSIFICATION <- factor(diversities$CLINICAL.CLASSIFICATION, 
                                levels = c("Control","Mild","Hospitalized"))

annotation <- data.frame(
  x = c(1,2,3),
  y = c(1,1,1),
  label = c(paste("n=",
                  nrow(diversities[diversities$CLINICAL.CLASSIFICATION=="Control",])),
            paste("n=",
                  nrow(diversities[diversities$CLINICAL.CLASSIFICATION=="Mild",])),
            paste("n=",
                  nrow(diversities[diversities$CLINICAL.CLASSIFICATION=="Hospitalized",])))
)


ggplot(diversities, aes(x = CLINICAL.CLASSIFICATION, y = shannon )) + 
  ylab("Shannon Diversity") +
  xlab("Disease Outcome") +
  labs(fill="Disease Outcome") +
  geom_boxplot(aes(fill=AGE_GROUP)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  geom_label(data=annotation, aes( x=x, y=y, label=label),
            color="black", 
            size=7 , angle=45, fontface="bold" )
  add_pvalue(df_p_val[df_p_val$p<0.06,], label = "p",
             step.increase = 0.03)
  add_pvalue(df_p_val[df_p_val$p.adj.signif!="ns",],
             label = "p.adj.signif",
             y.position = 2,
             tip.length = -0.02,
             step.increase = 0.03)


df_p_val <- rstatix::wilcox_test(diversities, simpson ~ CLINICAL.CLASSIFICATION) %>%
  rstatix::add_xy_position()


diversities$CLINICAL.CLASSIFICATION <- factor(diversities$CLINICAL.CLASSIFICATION, 
                                              levels = c("Control","Mild","Hospitalized"))

annotation <- data.frame(
  x = c(1,2,3),
  y = c(0.7,0.7,0.7),
  label = c(paste("n=",
                  nrow(diversities[diversities$CLINICAL.CLASSIFICATION=="Control",])),
            paste("n=",
                  nrow(diversities[diversities$CLINICAL.CLASSIFICATION=="Mild",])),
            paste("n=",
                  nrow(diversities[diversities$CLINICAL.CLASSIFICATION=="Hospitalized",])))
)


ggplot(diversities, aes(x = CLINICAL.CLASSIFICATION, y = simpson )) + 
  ylab("Simpson Diversity") +
  xlab("Disease Outcome") +
  labs(fill="Disease Outcome") +
  geom_boxplot(aes(fill=CLINICAL.CLASSIFICATION)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" )+
  add_pvalue(df_p_val[df_p_val$p<0.06,], label = "p",
             step.increase = 0.03)+
  add_pvalue(df_p_val[df_p_val$p.adj.signif!="ns",],
             label = "p.adj.signif",
             y.position = 0.75,
             tip.length = -0.02,
             step.increase = 0.03)


BH <- na.omit(diversities[(diversities$CITY=="BH") & (diversities$CLINICAL.CLASSIFICATION!="Moderate"),])

df_p_val <- rstatix::wilcox_test(BH, shannon ~ CLINICAL.CLASSIFICATION) %>%
  rstatix::add_xy_position()


BH$CLINICAL.CLASSIFICATION <- factor(BH$CLINICAL.CLASSIFICATION, 
                                              levels = c("Control","Mild","Severe"))

annotation <- data.frame(
  x = c(1,2,3),
  y = c(1,1,1),
  label = c(paste("n=",
                  nrow(BH[BH$CLINICAL.CLASSIFICATION=="Control",])),
            paste("n=",
                  nrow(BH[BH$CLINICAL.CLASSIFICATION=="Mild",])),
            paste("n=",
                  nrow(BH[BH$CLINICAL.CLASSIFICATION=="Severe",])))
)


ggplot(BH, aes(x = CLINICAL.CLASSIFICATION, y = shannon)) + 
  ylab("Shannon Diversity") +
  xlab("Disease Outcome") +
  labs(fill="Disease Outcome") +
  geom_boxplot(aes(fill=CLINICAL.CLASSIFICATION)) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"),
        legend.title=element_text(size=0),
        axis.title.x = element_text(size=0),
        legend.position = -1)
  add_pvalue(df_p_val[df_p_val$p<0.06,], label = "p",
             y.position = 7,
             step.increase = 0.05)
  add_pvalue(df_p_val[df_p_val$p.adj.signif!="ns",],
             label = "p.adj.signif",
             y.position = 2,
             tip.length = -0.02,
             step.increase = 0.03)


age_clin <- diversities %>%
  add_column(age_clin = NA)
age_clin$age_clin <- paste(age_clin$AGE_GROUP,age_clin$CLINICAL.CLASSIFICATION,sep=" ")

df_p_val <- rstatix::wilcox_test(age_clin, shannon ~ age_clin) %>%
  rstatix::add_xy_position()


ggplot(age_clin, aes(x = age_clin, y = shannon)) + 
  ylab("Shannon Diversity") +
  xlab("Disease - Age") +
  labs(fill="Disease Outcome") +
  geom_boxplot(aes(fill=AGE_GROUP)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" )+
  add_pvalue(df_p_val[df_p_val$p<0.2,], label = "p",
             step.increase = 0.03)
  add_pvalue(df_p_val[df_p_val$p.adj.signif!="ns",],
             label = "p.adj.signif",
             y.position = 2,
             tip.length = -0.02,
             step.increase = 0.03)


########## TESTE DE NORMALIDADE #############3
library("ggpubr")
ggdensity(diversities$shannon)
shapiro.test(diversities$shannon) #se p > 0.05 => normal 


teste <- diversities %>%
  add_column( OUTCOME_AGE=paste(diversities$CLINICAL.CLASSIFICATION,diversities$AGE_GROUP))


annotation <- data.frame() %>%
  add_column(x=NA) %>%
  add_column(y=NA) %>%
  add_column(label=NA) 
count =1
for(group in unique(teste$OUTCOME_AGE)){
  tmp <- data.frame(
    x=group,
    y=0.5,
    label=paste("n=",
                  nrow(teste[teste$OUTCOME_AGE==group,])))
  
  count = count+1
  annotation <- rbind(annotation,tmp)
}

df_p_val <- rstatix::wilcox_test(teste, simpson ~ OUTCOME_AGE) %>%
  rstatix::add_xy_position()


ggplot(teste, aes(x = OUTCOME_AGE, y = simpson)) + 
  ylab("simpson Diversity") +
  xlab("Disease - Age") +
  labs(fill="Age Group") +
  geom_boxplot(aes(fill=AGE_GROUP), width=0.5, position=position_dodge(0.5)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  ggtitle("Top100 diversity Age ˜ Disease outcome") +
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" ) +
  add_pvalue(df_p_val[df_p_val$p<0.05,], label = "p",
             step.increase = 0.04,
             y.position = 1) 


add_pvalue(df_p_val[df_p_val$p.adj.signif!="ns",],
           label = "p.adj.signif",
           y.position = 2,
           tip.length = -0.02,
           step.increase = 0.03)






teste <- diversities %>%
  add_column( OUTCOME_CITY=paste(diversities$CITY,diversities$CLINICAL.CLASSIFICATION))


annotation <- data.frame() %>%
  add_column(x=NA) %>%
  add_column(y=NA) %>%
  add_column(label=NA) 
count =1
for(group in unique(teste$OUTCOME_CITY)){
  tmp <- data.frame(
    x=group,
    y=1,
    label=paste("n=",
                nrow(teste[teste$OUTCOME_CITY==group,])))
  
  count = count+1
  annotation <- rbind(annotation,tmp)
}

df_p_val <- rstatix::wilcox_test(teste, shannon ~ OUTCOME_CITY) %>%
  rstatix::add_xy_position()


ggplot(teste, aes(x = OUTCOME_CITY, y = shannon)) + 
  ylab("Shannon Diversity") +
  xlab("Disease - Age") +
  labs(fill="City") +
  geom_boxplot(aes(fill=CITY), width=0.5, position=position_dodge(0.5)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  ggtitle("Top100 diversity City ˜ Disease outcome") +
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" )+
  add_pvalue(df_p_val[df_p_val$p<0.05,], label = "p",
             step.increase = 0.04,
             y.position = 5)  


add_pvalue(df_p_val[df_p_val$p.adj.signif!="ns",],
           label = "p.adj.signif",
           y.position = 2,
           tip.length = -0.02,
           step.increase = 0.03)


teste <- diversities[diversities$CITY=="BH" | diversities$CITY=="Control",] %>%
  add_column( OUTCOME_AGE=paste(diversities[diversities$CITY=="BH" | diversities$CITY=="Control",]$CLINICAL.CLASSIFICATION,diversities[diversities$CITY=="BH" | diversities$CITY=="Control",]$AGE_GROUP))


annotation <- data.frame() %>%
  add_column(x=NA) %>%
  add_column(y=NA) %>%
  add_column(label=NA) 
count =1
for(group in unique(teste$CLINICAL.CLASSIFICATION)){
  tmp <- data.frame(
    x=group,
    y=1,
    label=paste("n=",
                nrow(teste[teste$CLINICAL.CLASSIFICATION==group,])))
  
  count = count+1
  annotation <- rbind(annotation,tmp)
}

df_p_val <- rstatix::wilcox_test(teste, shannon ~ CLINICAL.CLASSIFICATION) %>%
  rstatix::add_xy_position()


ggplot(teste, aes(x = CLINICAL.CLASSIFICATION, y = shannon)) + 
  ylab("Shannon Diversity") +
  xlab("Disease Outcome") +
  labs(fill="Outcome") +
  geom_boxplot(aes(fill=CLINICAL.CLASSIFICATION), width=0.5, position=position_dodge(0.5)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) 
  add_pvalue(df_p_val[df_p_val$p<0.05,], label = "p",
             step.increase = 0.04,
             y.position = 5) 
  ggtitle("BH only TOP100 diversity Disease outcome") +
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" )
 


########## MILD ##########

mild <- na.omit(diversities[(((diversities$CITY=="BH") & (diversities$CLINICAL.CLASSIFICATION=="Mild")) | (diversities$CITY=="GV")) & (diversities$CLINICAL.CLASSIFICATION=="Mild"),])

teste <- mild %>%
  add_column( OUTCOME_AGE=paste(diversities[diversities$CITY=="BH" | diversities$CITY=="Control",]$CLINICAL.CLASSIFICATION,diversities[diversities$CITY=="BH" | diversities$CITY=="Control",]$AGE_GROUP))


annotation <- data.frame() %>%
  add_column(x=NA) %>%
  add_column(y=NA) %>%
  add_column(label=NA) 
count =1
for(group in unique(mild$CITY)){
  tmp <- data.frame(
    x=group,
    y=1,
    label=paste("n=",
                nrow(mild[mild$CITY==group,])))
  
  count = count+1
  annotation <- rbind(annotation,tmp)
}

df_p_val <- rstatix::t_test(mild, shannon ~ CITY) %>%
  rstatix::add_xy_position()

df_p_val <- df_p_val %>%
  add_column(p_ast="**")


ggplot(mild, aes(x = CITY, y = shannon)) + 
  ylab("Shannon Diversity") +
  labs(fill="Outcome") +
  geom_boxplot(aes(fill=CLINICAL.CLASSIFICATION), width=0.5, position=position_dodge(0.5)) +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"),
        legend.title=element_text(size=0),
        axis.title.x = element_text(size=0),
        legend.position = -1) +
  add_pvalue(df_p_val[df_p_val$p<0.05,], label = "p_ast",
             step.increase = 0.04,
             y.position = 7) 
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" )
  

teste <- diversities[(diversities$CITY=="BH") & (diversities$CLINICAL.CLASSIFICATION!="Moderate"),]


annotation <- data.frame() %>%
  add_column(x=NA) %>%
  add_column(y=NA) %>%
  add_column(label=NA) 
count =1
for(group in unique(teste$CLINICAL.CLASSIFICATION)){
  tmp <- data.frame(
    x=group,
    y=1,
    label=paste("n=",
                nrow(teste[teste$CLINICAL.CLASSIFICATION==group,])))
  
  count = count+1
  annotation <- rbind(annotation,tmp)
}

df_p_val <- rstatix::wilcox_test(teste, shannon ~ CLINICAL.CLASSIFICATION) %>%
  rstatix::add_xy_position()


ggplot(teste, aes(x = CLINICAL.CLASSIFICATION, y = shannon)) + 
  ylab("Shannon Diversity") +
  xlab("Disease Outcome") +
  labs(fill="Outcome") +
  geom_boxplot(aes(fill=CLINICAL.CLASSIFICATION), width=0.5, position=position_dodge(0.5)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.title=element_text(size=12)) +
  ggtitle("BH only D50 Disease outcome")+
  geom_label(data=annotation, aes( x=x, y=y, label=label),
             color="black", 
             size=7 , angle=45, fontface="bold" )
  add_pvalue(df_p_val[df_p_val$p<0.05,], label = "p",
             step.increase = 0.04,
             y.position = 5) 
