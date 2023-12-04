library(dplyr)
library(tidyr)
library(ggplot2)



if (!(dir.exists("../output/PatchSeqDataset/Electrophys_feature_plots"))){
  dir.create("../output/PatchSeqDataset/Electrophys_feature_plots")
}

if (!(dir.exists("../output/PatchClampDataset/Electrophys_feature_plots"))){
  dir.create("../output/PatchClampDataset/Electrophys_feature_plots")
}

########## PatchSeq Dataset ############

#load met metadata
met <- read.csv("../data/PatchSeqDataset/PatchSeq_metadata.csv")
t_labels <- met[,c("cell_specimen_id","transcriptomics_sample_id", "T.type.Label")]
t_labels$transcriptomics_sample_id <- gsub("-",".", t_labels$transcriptomics_sample_id)
t_labels_p <- separate(t_labels, col = T.type.Label, into = c("type" , "subtype1", "subtype2" ), sep = " ")
t_labels_p <- unite(t_labels_p, type, subtype1, col = "subtype_name", sep = " ", remove = FALSE)

#load ephys metadata
ephys <- read.csv("../data/PatchSeqDataset/PatchSeq_EP_features.csv")

met_short <- t_labels_p %>% select(cell_specimen_id, type,transcriptomics_sample_id) %>% arrange(cell_specimen_id) %>%  mutate(type = ifelse(type=="", "Unknown", type))
ephys_short <- ephys %>%  select(cell_specimen_id ,peak_t_ramp, threshold_t_ramp, threshold_i_ramp,threshold_v_ramp, f_i_curve_slope, adaptation, avg_isi) %>%  arrange(cell_specimen_id) %>%  mutate(AP_half_width = (peak_t_ramp - threshold_t_ramp)*10^3/2)

#loading K=3 NSS cluster results
clus <- read.csv("../output/PatchSeqDataset/NSS_clusters/NSS_clusters_k3.csv", header = FALSE)
clus$V1 <- gsub("-",".", clus$V1)
clus <- clus %>%  arrange(V1)
colnames(clus)[1] <- "transcriptomics_sample_id"

all_data <- left_join(ephys_short, met_short, by = "cell_specimen_id")
all_data <- all_data[all_data$transcriptomics_sample_id %in% clus$transcriptomics_sample_id,]
all_data <- all_data %>%  arrange(transcriptomics_sample_id)
all_data <- left_join(all_data,clus, by = "transcriptomics_sample_id")
all_data <- all_data %>% mutate(V2 = as.factor(V2))




#all_data$cell_specimen_id <- factor(all_data$cell_specimen_id, levels = all_data$cell_specimen_id)
all_data <- all_data %>%  mutate(type = ifelse(type %in%  c("Serpinf1", "Sncg", "Meis2"),"Vip", type)) %>%  arrange(type, AP_half_width)
levels(all_data$V2) <- c("EC1","EC0", "EC2")
all_data$V2 <-factor(all_data$V2, levels= c("EC0","EC1", "EC2"))


all_data %>%  ggplot(aes(x=V2, y =AP_half_width/1000)) + geom_violin(aes(fill=V2)) + scale_fill_manual(values = c("EC0" = "#440154", "EC1" = "#21908c", "EC2" = "#fde725"))+ 
  geom_boxplot( alpha = 0.5) + ggtitle(expression("AP"[italic(halfwidth)]*" (ms)"))+
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 55) ,axis.title = element_blank(), 
        plot.title = element_text(size = 65, hjust = 0.5, face = "bold")) + scale_y_continuous(limits = c(0, 1.01))
ggsave(path = "../output/PatchSeqDataset/Electrophys_feature_plots", filename = "AP_Patch_seq.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)


all_data %>%  ggplot(aes(x=V2, y =f_i_curve_slope)) + geom_violin(aes(fill=V2)) + scale_fill_manual(values = c("EC0" = "#440154", "EC1" = "#21908c", "EC2" = "#fde725"))+ 
  geom_boxplot( alpha = 0.5) + ggtitle("f-I curve slope")+
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 55) ,axis.title = element_blank(), 
        plot.title = element_text(size = 65, hjust = 0.5)) 
ggsave(path = "../output/PatchSeqDataset/Electrophys_feature_plots", filename = "ficurve_Patch_seq.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

all_data %>%  ggplot(aes(x=V2, y =avg_isi)) + geom_violin(aes(fill=V2)) + scale_fill_manual(values = c("EC0" = "#440154", "EC1" = "#21908c", "EC2" = "#fde725"))+ 
  geom_boxplot( alpha = 0.5) + ggtitle("Average ISI (ms)")+
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 55) ,axis.title = element_blank(), 
        plot.title = element_text(size = 65, hjust = 0.5)) 
ggsave(path = "../output/PatchSeqDataset/Electrophys_feature_plots", filename = "avg_ISI_Patch_seq.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

all_data %>%  ggplot(aes(x=V2, y =threshold_i_ramp)) + geom_violin(aes(fill=V2)) + scale_fill_manual(values = c("EC0" = "#440154", "EC1" = "#21908c", "EC2" = "#fde725"))+ 
  geom_boxplot( alpha = 0.5) + ggtitle(expression("I"[italic(THR)]*" (pA)"))+
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 55) ,axis.title = element_blank(), 
        plot.title = element_text(size = 65, hjust = 0.5, face = "bold")) 
ggsave(path = "../output/PatchSeqDataset/Electrophys_feature_plots", filename = "th_I_clus_Patch_seq.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

all_data %>%  ggplot(aes(x=type, y =threshold_i_ramp)) + geom_violin(aes(fill=type)) + scale_fill_manual(values = c("Pvalb" = "#53B400", "Sst" = "#A58AFF", "Vip" = "#FB61D7", "Lamp5" = "#F8766D"))+ 
  #geom_boxplot( alpha = 0.5) +
  ggtitle(expression("I"[italic(THR)]*" (pA)"))+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 30) ,axis.title = element_blank(), 
        plot.title = element_text(size = 40, hjust = 0.5, face = "bold"),legend.text = element_text(size=40), legend.title = element_text(size = 40) ) +
  geom_jitter(aes(color=V2), size = 2) + scale_color_manual(values = c("EC0" = "#440154", "EC1" = "#21908c", "EC2" = "#fde725"))+ guides(fill=FALSE) + labs(color = "Cluster")
ggsave(path = "../output/PatchSeqDataset/Electrophys_feature_plots", filename = "th_I_type_puntini_Patch_seq.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)



medians <- all_data %>% group_by(V2) %>%  summarise(AP_median = median(AP_half_width), fi_median = median(f_i_curve_slope), ISI_median = median(avg_isi), T_v = median(threshold_v_ramp), T_I= median(threshold_i_ramp))
write.csv(medians, "../output/PatchSeqDataset/Electrophys_feature_plots/medians.csv")


################ PatchClamp Dataset ##############

#load cre metadata
met_pc <- read.csv("../data/PatchClampDataset/PatchClamp_metadata.csv", sep = ";")
t_labels <- met_pc[,c("cell_specimen_id","transcriptomics_sample_id", "T.type.Label")]
t_labels$transcriptomics_sample_id <- gsub("-",".", t_labels$transcriptomics_sample_id)
t_labels_p <- separate(t_labels, col = T.type.Label, into = c("type" , "subtype1", "subtype2" ), sep = " ")
t_labels_p <- unite(t_labels_p, type, subtype1, col = "subtype_name", sep = " ", remove = FALSE)

#load ephys metadata
ephys_pc <- read.csv("../data/PatchClampDataset/PatchClamp_EP_features.csv", sep = ",")

clus_pc_p <- read.csv("../output/PatchClampDataset/NSS_clusters/NSS_clusters_k2.csv", header = FALSE)
line_0_p <- read.csv("../output/PatchClampDataset/cre_lines_clusters/lines_k2_cl0.csv", header = FALSE)
line_1_p <- read.csv("../output/PatchClampDataset/cre_lines_clusters/lines_k2_cl1.csv", header = FALSE)

line_all_p <- rbind(line_0_p, line_1_p) 
line_clus_p <- cbind(line_all_p, clus_pc_p)
line_clus_p$V1 <- NULL
colnames(line_clus_p)<- c("specimen_id", "line", "cluster")
line_clus_p <- line_clus_p %>% arrange(specimen_id)
type_pc <- read.csv("../output/PatchClampDataset/MISC/cre_lines_type.csv", header = FALSE) %>%  arrange(V2)
colnames(type_pc)<- c("type","specimen_id")
type_pc <- left_join(line_clus_p,type_pc, by="specimen_id")
ephys_pc_p <- ephys_pc %>%  filter(specimen_id %in% type_pc$specimen_id ) %>%  select(specimen_id , peak_t_ramp, threshold_t_ramp, threshold_i_ramp,threshold_v_ramp, f_i_curve_slope, avg_isi) %>% mutate(AP_half_width = (peak_t_ramp - threshold_t_ramp)*10^6/2) %>%   arrange(specimen_id)




PC_DATA_2 <- left_join(ephys_pc_p, type_pc, by= "specimen_id") %>% mutate(cluster = as.factor(cluster))
levels(PC_DATA_2$cluster) <- c("EC1", "EC0")
PC_DATA_2$cluster <-factor(PC_DATA_2$cluster, levels= c("EC0","EC1"))



PC_DATA_2 %>%  ggplot(aes(x=cluster, y =AP_half_width/1000)) + geom_violin(aes(fill=cluster)) + scale_fill_manual(values = c("EC0" = "#440154", "EC2" = "#21908c", "EC1" = "#fde725"))+ 
  geom_boxplot( alpha = 0.5) + ggtitle(expression("AP"[italic(halfwidth)]*" (ms)"))+
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 55) ,axis.title = element_blank(), 
        plot.title = element_text(size = 65, hjust = 0.5, face = "bold")) + scale_y_continuous(limits = c(0, 1.01))
ggsave(path = "../output/PatchClampDataset/Electrophys_feature_plots", filename = "AP_Patch_clamp.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

PC_DATA_2 %>%  ggplot(aes(x=cluster, y =f_i_curve_slope)) + geom_violin(aes(fill=cluster)) + scale_fill_manual(values = c("EC0" = "#440154", "EC2" = "#21908c", "EC1" = "#fde725"))+ 
  geom_boxplot( alpha = 0.5) + ggtitle("f-I curve slope")+
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 55) ,axis.title = element_blank(), 
        plot.title = element_text(size = 65, hjust = 0.5)) 
ggsave(path = "../output/PatchClampDataset/Electrophys_feature_plots", filename = "ficurve_Patch_clamp.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

PC_DATA_2 %>%  ggplot(aes(x=cluster, y =avg_isi)) + geom_violin(aes(fill=cluster)) + scale_fill_manual(values = c("EC0" = "#440154", "EC2" = "#21908c", "EC1" = "#fde725"))+ 
  geom_boxplot( alpha = 0.5) + ggtitle("Average ISI (ms)")+
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 55) ,axis.title = element_blank(), 
        plot.title = element_text(size = 65, hjust = 0.5)) 
ggsave(path = "../output/PatchClampDataset/Electrophys_feature_plots", filename = "avg_ISI_Patch_clamp.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

PC_DATA_2 %>%  ggplot(aes(x=cluster, y =threshold_i_ramp)) + geom_violin(aes(fill=cluster)) + scale_fill_manual(values = c("EC0" = "#440154", "EC2" = "#21908c", "EC1" = "#fde725"))+ 
  geom_boxplot( alpha = 0.5) + ggtitle(expression("I"[italic(THR)]*" (pA)"))+
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 55) ,axis.title = element_blank(), 
        plot.title = element_text(size = 65, hjust = 0.5, face = "bold")) 
ggsave(path = "../output/PatchClampDataset/Electrophys_feature_plots", filename = "th_I_clus_Patch_clamp.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

PC_DATA_2 %>%  ggplot(aes(x=cluster, y =threshold_v_ramp)) + geom_violin(aes(fill=cluster)) + scale_fill_manual(values = c("EC0" = "#440154", "EC1" = "#21908c", "EC2" = "#fde725"))+ 
  geom_boxplot( alpha = 0.5) + ggtitle("Voltage Threshold (pA)")+
  theme(legend.position = "none", axis.title.x=element_blank(),
        axis.ticks.x=element_blank(), axis.text = element_text(size = 30) ,axis.title = element_blank(), 
        plot.title = element_text(size = 40, hjust = 0.5, face = "bold")) 
ggsave(path = "../output/PatchClampDataset/Electrophys_feature_plots", filename = "th_V_clus_Patch_clamp.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

medians_pc <- PC_DATA_2 %>% group_by(cluster) %>%  summarise(AP_median = median(AP_half_width), fi_median = median(f_i_curve_slope), ISI_median = median(avg_isi), T_v = median(threshold_v_ramp), T_I= median(threshold_i_ramp))
write.csv(medians_pc, "../output/PatchClampDataset/Electrophys_feature_plots/medians_pc.csv")

