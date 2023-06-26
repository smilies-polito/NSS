#
# Copyright © 2022 Politecnico di Torino, Control and Computer Engineering Department, SMILIES group
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
# associated documentation files (the "Software"), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial 
# portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT 
# LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

#######################################################################################
#
# LOAD REQUIRED PACKAGES
#
#######################################################################################

library(data.table)
library(readsparse)
library(aricode)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(stringr)
library(patchwork)

#######

if (!(dir.exists("../output"))){
  dir.create("../output")
}
if (!(dir.exists("../output/PatchSeqDataset"))){
  dir.create("../output/PatchSeqDataset")
}
if (!(dir.exists("../output/PatchSeqDataset/IMAGES"))){
    dir.create("../output/PatchSeqDataset/IMAGES")
}
if (!(dir.exists("../output/PatchSeqDataset/GO_enrichment_analysis"))){
    dir.create("../output/PatchSeqDataset/GO_enrichment_analysis")
}
#########

args <- commandArgs(trailingOnly = TRUE)

# Check if argument was passed
if (length(args) > 0) {
    K <- args[1]
    print(paste("Processing for K =", K))
} else {
    print("No arguments provided.")
}

# Load the data and preprocess
print("Data Loading")
Patch_seq <- read.csv("../data/PatchSeqDataset/count.csv", header=TRUE, row.names=1)
print("Data processing")
Patch_seq_S <- CreateSeuratObject(Patch_seq)
Patch_seq_S <- NormalizeData(Patch_seq_S)
Patch_seq_S <- FindVariableFeatures(Patch_seq_S)
Patch_seq_S <- ScaleData(Patch_seq_S)
Patch_seq_S <- RunPCA(Patch_seq_S)
Patch_seq_S <- FindNeighbors(Patch_seq_S)
Patch_seq_S <- FindClusters(Patch_seq_S, resolution = 0.8)
Patch_seq_S <- RunUMAP(Patch_seq_S, dims = 1:10)

# Adding metadata and subsetting to inhibitory neurons
met_labels <- read.csv("../data/PatchSeqDataset/PatchSeq_metadata.csv")
t_labels <- met_labels[,c("transcriptomics_sample_id", "T.type.Label")]
t_labels$transcriptomics_sample_id <- gsub("-",".", t_labels$transcriptomics_sample_id)
t_labels <- t_labels[t_labels$transcriptomics_sample_id %in% colnames(Patch_seq_S),]
t_labels <- t_labels[match(colnames(Patch_seq_S), t_labels$transcriptomics_sample_id),]
t_labels_p <- separate(t_labels, col = T.type.Label, into = c("type" , "subtype1", "subtype2" ), sep = " ")
t_labels_p <- unite(t_labels_p, type, subtype1, sep = " ",col = "subtype_name", remove = FALSE)

#Patch_seq_sub <- subset(Patch_seq_S, cells = t_labels$transcriptomics_sample_id)

met_labels <- met_labels[,c("transcriptomics_sample_id", "MET.type.Label")]
met_labels$transcriptomics_sample_id <- gsub("-",".", met_labels$transcriptomics_sample_id)
met_labels <- met_labels[met_labels$transcriptomics_sample_id %in% colnames(Patch_seq_S),]
met_labels <- met_labels[match(colnames(Patch_seq_S), met_labels$transcriptomics_sample_id),]


Patch_seq_S <- AddMetaData(Patch_seq_S, t_labels$T.type.Label, col.name = "T_labels")
Patch_seq_S <- AddMetaData(Patch_seq_S, t_labels_p$type, col.name = "T_labels_main")
Patch_seq_S <- AddMetaData(Patch_seq_S, met_labels$MET.type.Label, col.name = "MET_labels")

# Adding Eph clusters labels (both k=2 and k=3)
if (K == 3){
    fs_label <- read.csv("../output/PatchSeqDataset/NSS_clusters/NSS_clusters_k3.csv", header=FALSE)
    fs_label$V1 <- gsub("-", ".", fs_label$V1)
    fs_cl <- fs_label[fs_label$V1 %in% colnames(Patch_seq_S),]
    fs_cl$V2 <- paste0("EC",fs_cl$V2 )
    Patch_seq_sub_fs <- subset(Patch_seq_S, cells = fs_label$V1)
    Patch_seq_sub_fs <- AddMetaData(Patch_seq_sub_fs, fs_cl$V2, col.name = "Cluster_fs")
    
    Idents(Patch_seq_sub_fs) <- "Cluster_fs"
    new.cluster.ids <- c("EC2", "EC0", "EC1")
    names(new.cluster.ids) <- levels(Patch_seq_sub_fs)
    new.cluster.ids <- sort(new.cluster.ids)
    Patch_seq_sub_fs <- RenameIdents(Patch_seq_sub_fs ,new.cluster.ids)
    DimPlot(Patch_seq_sub_fs, pt.size =2.5, cols = c("#440154","#21908c","#fde725"))
    
    # Uncommnet this part to generate the plot with only one cluster highlighted as the final figure in the article
    
    #DimPlot(Patch_seq_sub_fs, pt.size =2, cols = c("#440154","#D3D3D3","#D3D3D3")) + NoLegend()
    #ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "purple_high.png", width = 640, height = 480, units= "px",scale = 3.5)
    #DimPlot(Patch_seq_sub_fs, pt.size =2, cols = c("#D3D3D3","#21908c","#D3D3D3")) + NoLegend()
    #ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "green_high.png", width = 640, height = 480, units= "px",scale = 3.5)
    #DimPlot(Patch_seq_sub_fs, pt.size = 2, cols = c("#D3D3D3","#D3D3D3","#fde725")) + NoLegend()
    #ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "yellow_high.png", width = 640, height = 480, units= "px",scale = 3.5)
    
    
    DimPlot(Patch_seq_sub_fs, pt.size =2, cols = c("#440154","#21908c","#fde725")) + ggtitle("Eph Clusters (K=3)") +
        theme(legend.text = element_text(size = 50), axis.title = element_text(size=40), plot.title = element_text(size = 40, hjust = 0.5, face = "bold")) +
        guides(colour = guide_legend(override.aes = list(size=8))) 
    ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_embedding_NSS_labels_k3.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)
}


if (K==2){
    fs_label <- read.csv("../output/PatchSeqDataset/NSS_clusters/NSS_clusters_k2.csv", header=FALSE)
    fs_label$V1 <- gsub("-", ".", fs_label$V1)
    fs_cl_2 <- fs_label[fs_label$V1 %in% colnames(Patch_seq_S),]
    fs_cl_2$V2 <- paste0("EC",fs_cl_2$V2 )
    Patch_seq_sub_fs <- subset(Patch_seq_S, cells = fs_label$V1)
    Patch_seq_sub_fs <- AddMetaData(Patch_seq_sub_fs, fs_cl_2$V2, col.name = "Cluster_fs")
    
    Idents(Patch_seq_sub_fs) <- "Cluster_fs"
    new.cluster.ids <- c("EC1", "EC0")
    names(new.cluster.ids) <- levels(Patch_seq_sub_fs)
    new.cluster.ids <- sort(new.cluster.ids)
    Patch_seq_sub_fs <- RenameIdents(Patch_seq_sub_fs ,new.cluster.ids)
    DimPlot(Patch_seq_sub_fs, pt.size =2, cols = c( "#440154","#fde725")) + ggtitle("Eph Clusters (K=2)") +
        theme(legend.text = element_text(size = 50), axis.title = element_text(size=40), plot.title = element_text(size = 40, hjust = 0.5, face = "bold")) +
        guides(colour = guide_legend(override.aes = list(size=8))) 
    ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_embedding_NSS_labels_k2.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)
}

Idents(Patch_seq_sub_fs) <- "T_labels_main"

# for simplifying the visualization, "Serpinf1", "Sncg" and "Meis2" were labeled as "Vip"
# uncomment the following line to plot higher-resolution cell type labels:
# new.cluster.ids <- c("Lamp5", "Sst", "Vip", "Pvalb", "Serpinf1", "Sncg", "Meis2")
new.cluster.ids <- c("Lamp5", "Sst", "Vip", "Pvalb", "Vip", "Vip", "Vip")

names(new.cluster.ids) <- levels(Patch_seq_sub_fs)
new.cluster.ids <- sort(new.cluster.ids)
Patch_seq_sub_fs <- RenameIdents(Patch_seq_sub_fs, new.cluster.ids)

# add these colors to the cols parameter to support higher-resolution cell type labels including "Serpinf1", "Sncg" and "Meis2" "#FF5722", "#212121", "#AAAAAA"
DimPlot(Patch_seq_sub_fs, pt.size =2, cols = c("#F8766D","#53B400","#A58AFF","#FB61D7")) + ggtitle("Transcriptional metadata labels") +
    theme(legend.text = element_text(size = 50), axis.title = element_text(size=40), plot.title = element_text(size = 40, hjust = 0.5, face = "bold")) +
    guides(colour = guide_legend(override.aes = list(size=8)))
ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_embedding_transcriptomic_labels.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

g_p <- FeaturePlot(Patch_seq_sub_fs, features = c("Pvalb"), pt.size = 2 ) + theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text = element_text(size=20), axis.title = element_text(size=30)) 
g_s <- FeaturePlot(Patch_seq_sub_fs, features = c("Sst"), pt.size = 2 ) + theme(plot.title = element_text(hjust = 0.5, size = 30) , axis.text = element_text(size=20), axis.title = element_text(size=30))
g_v <- FeaturePlot(Patch_seq_sub_fs, features = c("Vip"), pt.size = 2 ) + theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text = element_text(size=20), axis.title = element_text(size=30))
g_l <- FeaturePlot(Patch_seq_sub_fs, features = c("Lamp5"), pt.size = 2 ) + theme(plot.title = element_text(hjust = 0.5, size = 30), axis.text = element_text(size=20), axis.title = element_text(size=30))
(g_p/g_v)+ (g_s/g_v)
wrap_plots(g_p,g_s,g_v,g_l, ncol = 2, nrow = 2)
ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_embedding_marker_expression.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)


########## Sst zoom ###########
print("Focus on Sst")
Idents(Patch_seq_sub_fs) <- "T_labels_main"
Patch_seq_Sst <- subset(Patch_seq_sub_fs, idents =  "Sst")
DimPlot(Patch_seq_Sst)
Patch_seq_Sst <- NormalizeData(Patch_seq_Sst)
Patch_seq_Sst <- FindVariableFeatures(Patch_seq_Sst)
Patch_seq_Sst <- ScaleData(Patch_seq_Sst)
Patch_seq_Sst <- RunPCA(Patch_seq_Sst)
Patch_seq_Sst <- FindNeighbors(Patch_seq_Sst)
Patch_seq_Sst <- FindClusters(Patch_seq_Sst, resolution = 0.5)
Patch_seq_Sst <- RunUMAP(Patch_seq_Sst, 1:20)
#colnames(Patch_seq_Sst)
DimPlot(Patch_seq_Sst, group.by = "T_labels", pt.size = 2, )+ ggtitle("Sst subtype labels") + theme(legend.text = element_text(size = 35), axis.title = element_text(size=40), title = element_text(size = 40)) +
    guides(colour = guide_legend(override.aes = list(size=8), ncol = 1)) 

ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_embedding_transcriptomic_labels_Sst_subset.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)

if(K==2){
    Idents(Patch_seq_Sst) <- "Cluster_fs"
    new.cluster.ids <- c("EC1", "EC0")
    names(new.cluster.ids) <- levels(Patch_seq_Sst)
    new.cluster.ids <- sort(new.cluster.ids)
    Patch_seq_Sst <- RenameIdents(Patch_seq_Sst ,new.cluster.ids)
    DimPlot(Patch_seq_Sst, pt.size = 2, cols = c("#440154","#fde725") ) + ggtitle("Eph Clusters (K=3)") +
        theme(legend.text = element_text(size = 50), axis.title = element_text(size=40), plot.title = element_text(size = 40, hjust = 0.5, face = "bold")) +
        guides(colour = guide_legend(override.aes = list(size=8))) 
    ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_embedding_NSS_labels_k2_Sst_subset.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)
}
if(K==3){
    Idents(Patch_seq_Sst) <- "Cluster_fs"
    new.cluster.ids <- c("EC2", "EC0", "EC1")
    names(new.cluster.ids) <- levels(Patch_seq_Sst)
    new.cluster.ids <- sort(new.cluster.ids)
    Patch_seq_Sst <- RenameIdents(Patch_seq_Sst ,new.cluster.ids)
    #Patch_seq_Sst <- RenameIdents(Patch_seq_Sst ,new.cluster.ids)
    DimPlot(Patch_seq_Sst, pt.size = 2, cols = c("#440154","#21908c","#fde725") ) + ggtitle("Eph Clusters (K=3)") +
        theme(legend.text = element_text(size = 50), axis.title = element_text(size=40), plot.title = element_text(size = 40, hjust = 0.5, face = "bold")) +
        guides(colour = guide_legend(override.aes = list(size=8))) 
    ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_embedding_NSS_labels_k3_Sst_subset.pdf", width = 1920, height = 1080, units= "px",scale = 3.5)
}

#################

# Differential expression
print("DE analysis")
if (K==3){
    Idents(Patch_seq_sub_fs) <- "Cluster_fs"
    new.cluster.ids <- c("EC2", "EC0", "EC1")
    names(new.cluster.ids) <- levels(Patch_seq_sub_fs)
    new.cluster.ids <- sort(new.cluster.ids)
    Patch_seq_sub_fs <- RenameIdents(Patch_seq_sub_fs ,new.cluster.ids)
    
    EC0_vs_EC1 <- FindMarkers(Patch_seq_sub_fs, ident.1 = "EC0", ident.2 = "EC1")
    EC0_vs_EC1$gene <- rownames(EC0_vs_EC1)
    EC0_vs_EC1_P <- EC0_vs_EC1 %>% filter(avg_log2FC > 0)
    write.csv(EC0_vs_EC1_P, file = paste0("../output/PatchSeqDataset/GO_enrichment_analysis/NSS_cluster0_k",K,"_DE_genes.csv"))
    
    EC1_vs_EC2 <- FindMarkers(Patch_seq_sub_fs, ident.1 = "EC1", ident.2 = "EC2")
    EC1_vs_EC2$gene <- rownames(EC1_vs_EC2)
    EC1_vs_EC2_P <- EC1_vs_EC2 %>% filter(avg_log2FC > 0)
    
    EC0_vs_EC2 <- FindMarkers(Patch_seq_sub_fs, ident.1 = "EC0", ident.2 = "EC2")
    EC0_vs_EC2$gene <- rownames(EC0_vs_EC2)
    gene_01 <- EC0_vs_EC1 %>% filter(avg_log2FC > 0) %>% filter(str_detect(.$gene, "Kc")) %>% .$gene
    gene_12 <- EC1_vs_EC2 %>% filter(avg_log2FC > 0) %>% filter(str_detect(.$gene, "Kc")) %>% .$gene
    gene_pos <- gene_01[gene_01 %in% gene_12]
    
    gene_01 <- EC0_vs_EC1 %>% filter(avg_log2FC < 0) %>% filter(str_detect(.$gene, "Kc")) %>% .$gene
    gene_12 <- EC1_vs_EC2 %>% filter(avg_log2FC < 0) %>% filter(str_detect(.$gene, "Kc")) %>% .$gene
    gene_neg <- gene_01[gene_01 %in% gene_12]
    
    VlnPlot(Patch_seq_sub_fs, features = c("Kcnc2"), cols = c("#440154","#21908c","#fde725")) + NoLegend() + theme(axis.title.x=element_blank(),
                                                                                                                   axis.text.x=element_blank(),
                                                                                                                   axis.ticks.x=element_blank(), plot.title = element_blank())
    ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_Kcnc2_violin_plot_k3.pdf", width = 960, height = 260, units= "px",scale = 3.5)
    
    VlnPlot(Patch_seq_sub_fs, features = c("Kcnn2"), cols = c("#440154","#21908c","#fde725")) + NoLegend() + theme(axis.title.x=element_blank(),
                                                                                                                   axis.text.x=element_blank(),
                                                                                                                   axis.ticks.x=element_blank(), plot.title = element_blank())
    ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_Kcnn2_violin_plot_k3.pdf", width = 960, height = 260, units= "px",scale = 3.5)
    
}

if(K==2){
    Idents(Patch_seq_sub_fs) <- "Cluster_fs"
    new.cluster.ids <- c("EC1", "EC0")
    names(new.cluster.ids) <- levels(Patch_seq_sub_fs)
    new.cluster.ids <- sort(new.cluster.ids)
    Patch_seq_sub_fs <- RenameIdents(Patch_seq_sub_fs ,new.cluster.ids)
    EC0_vs_EC1 <- FindMarkers(Patch_seq_sub_fs, ident.1 = "EC0", ident.2 = "EC1")
    EC0_vs_EC1$gene <- rownames(EC0_vs_EC1)
    EC0_vs_EC1_P <- EC0_vs_EC1 %>% filter(avg_log2FC > 0)
    write.csv(EC0_vs_EC1_P, file = paste0("../output/PatchSeqDataset/GO_enrichment_analysis/NSS_cluster0_k",K,"_DE_genes.csv"))
    VlnPlot(Patch_seq_sub_fs, features = c("Kcnc2"), cols = c("#440154","#fde725")) + NoLegend() + theme(axis.title.x=element_blank(),
                                                                                                                   axis.text.x=element_blank(),
                                                                                                                   axis.ticks.x=element_blank(), plot.title = element_blank())
    ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_Kcnc2_violin_plot_k2.pdf", width = 960, height = 260, units= "px",scale = 3.5)
    
    VlnPlot(Patch_seq_sub_fs, features = c("Kcnn2"), cols = c("#440154","#fde725")) + NoLegend() + theme(axis.title.x=element_blank(),
                                                                                                                   axis.text.x=element_blank(),
                                                                                                                   axis.ticks.x=element_blank(), plot.title = element_blank())
    ggsave(path = "../output/PatchSeqDataset/IMAGES/", filename = "transcriptomic_Kcnn2_violin_plot_k2.pdf", width = 960, height = 260, units= "px",scale = 3.5)
    
}
############
