library(Seurat)
library(patchwork)
library(sctransform)
library(glmGamPoi)
library(future)
library(dplyr)
library(ggplot2)

setwd("/opt/home/buckcenter.org/fwu/scRNAseq_analysis/Sex_cellTypes_stim")

# change the current plan to access parallelization
plan("multisession", workers = 4)
options(future.globals.maxSize= +Inf) # increase maxSize of future function

Idents(Sample.combined)<-"orig.ident"
# # gene DEGs/markers of uG vs 1G stimulated PBMCs in Female
# markers.uGvs1G.stim.F <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="F.uG_stimulated", ident.2 ="F.1G_stimulated",
#                                min.pct=0.005, logfc.threshold = 0.1,
#                                test.use = "MAST")
# markers.uGvs1G.stim.F<-markers.uGvs1G.stim.F[order(markers.uGvs1G.stim.F$avg_log2FC, decreasing = T),] # sort by log2FC

# get a table for all the genes in uG vs 1G stimulated PBMCs in Female
markers.uGvs1G.stim.all.F <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="F.uG_stimulated", ident.2 ="F.1G_stimulated",
                                           min.pct=0, logfc.threshold = 0,
                                           test.use = "MAST")
markers.uGvs1G.stim.all.F<-markers.uGvs1G.stim.all.F[order(markers.uGvs1G.stim.all.F$avg_log2FC, decreasing = T),] # sort by log2FC

# save the marker lists of Female
# write.csv(markers.uGvs1G.stim.F, "Female_markers_uGvs1G_stimulated.csv")
write.csv(markers.uGvs1G.stim.all.F, "Overall_uGvs1G_Female.csv")

# get a table for all the genes in uG vs 1G stimulated PBMCs in Male
markers.uGvs1G.stim.all.M <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="M.uG_stimulated", ident.2 ="M.1G_stimulated",
                                           min.pct=0, logfc.threshold = 0,
                                           test.use = "MAST")
markers.uGvs1G.stim.all.M<-markers.uGvs1G.stim.all.M[order(markers.uGvs1G.stim.all.M$avg_log2FC, decreasing = T),] # sort by log2FC

# save the marker lists of Female
# write.csv(markers.uGvs1G.stim.M, "Female_markers_uGvs1G_stimulated.csv")
write.csv(markers.uGvs1G.stim.all.M, "Overall_uGvs1G_Male.csv")


# DEGs for each cell types
#Sample.combined$cell_sex <- paste0(Sample.combined$predicted.celltype.l2,"_",Sample.combined$orig.ident)
Idents(Sample.combined) <- "cell_sex"
cell4DEG<-unique(Sample.combined.stim$predicted.celltype.l2)
# for (d in cell4DEG[!cell4DEG %in% c("Unknow","ASDC","B memory","CD16 Mono","CD4 TEM","cDC1","CD4 Proliferating","ILC","NK Proliferating","Plasmablast")]){
#   DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_F.uG_stimulated"), ident.2 =paste0(d,"_F.1G_stimulated"))
#   write.csv(DEG.uG, paste0(d,"_uGvs1G_Female.csv"))
#   DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_M.uG_stimulated"), ident.2 =paste0(d,"_M.1G_stimulated"))
#   write.csv(DEG.uG, paste0(d,"_uGvs1G_Male.csv"))
# }

# full gene list
for (d in cell4DEG[!cell4DEG %in% c("Unknow","ASDC","B memory","CD16 Mono","CD4 TEM","cDC1","CD4 Proliferating","ILC","NK Proliferating","Plasmablast","CD4 Naive","pDC")]){
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_F.uG_stimulated"), ident.2 =paste0(d,"_F.1G_stimulated"),
                        min.pct=0, logfc.threshold = 0,
                        test.use = "MAST")
  write.csv(DEG.uG, paste0(d,"_uGvs1G_Female.csv"))
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_M.uG_stimulated"), ident.2 =paste0(d,"_M.1G_stimulated"),
                        min.pct=0, logfc.threshold = 0,
                        test.use = "MAST")
  write.csv(DEG.uG, paste0(d,"_uGvs1G_Male.csv"))
}

# compare the Female DEGs/Male DEGs
#setwd("../Sex_cellTypes_stim")
celltype4sex.stim <- gsub("_uGvs1G_Female.csv","",list.files(pattern="_Female.csv")) # grab the exist cell type

sex_DEGs.df <- data.frame()
for (d in celltype4sex.stim[!celltype4sex.stim %in% c("CD8 Naive","CD8 TCM","MAIT","NK_CD56bright")]){
  markers4sex.F <- read.csv(paste0(d,"_uGvs1G_Female.csv"))
  markers4sex.M <- read.csv(paste0(d,"_uGvs1G_Male.csv"))
  if (nrow(markers4dot) != 0){
    nDEGs_up.F <- nrow(filter(markers4sex.F,p_val_adj<0.05 & avg_log2FC>=0.1))
    nDEGs_down.F <- nrow(filter(markers4sex.F,p_val_adj<0.05 & avg_log2FC<=-0.1))
    nDEGs_up.M <- nrow(filter(markers4sex.M,p_val_adj<0.05 & avg_log2FC>=0.1))
    nDEGs_down.M <- nrow(filter(markers4sex.M,p_val_adj<0.05 & avg_log2FC<=-0.1))
    
    sex_DEGs.df.up <- data.frame(DEGs="Up-regulated",Cell_type=d,log2diffnDEGs=log2(nDEGs_up.F/nDEGs_up.M))
    sex_DEGs.df.down <- data.frame(DEGs="Down-regulated",Cell_type=d,log2diffnDEGs=log2(nDEGs_down.F/nDEGs_down.M))
    sex_DEGs.df <- rbind(sex_DEGs.df,sex_DEGs.df.up,sex_DEGs.df.down)
  }
}

celltype4sex.stim.other <- c("CD8 Naive","CD8 TCM","MAIT","NK_CD56bright")
nDEGs_up.F.sum <- 0
nDEGs_up.M.sum <- 0
nDEGs_down.F.sum <- 0
nDEGs_down.M.sum <- 0
for (o in celltype4sex.stim.other){
  markers4sex.F <- read.csv(paste0(o,"_uGvs1G_Female.csv"))
  markers4sex.M <- read.csv(paste0(o,"_uGvs1G_Male.csv"))
  
  nDEGs_up.F <- nrow(filter(markers4sex.F,p_val_adj<0.05 & avg_log2FC>=0.1))
  nDEGs_down.F <- nrow(filter(markers4sex.F,p_val_adj<0.05 & avg_log2FC<=-0.1))
  nDEGs_up.M <- nrow(filter(markers4sex.M,p_val_adj<0.05 & avg_log2FC>=0.1))
  nDEGs_down.M <- nrow(filter(markers4sex.M,p_val_adj<0.05 & avg_log2FC<=-0.1))
  
  nDEGs_up.F.sum <- nDEGs_up.F.sum + nDEGs_up.F
  nDEGs_up.M.sum <- nDEGs_up.M.sum + nDEGs_up.M
  nDEGs_down.F.sum <- nDEGs_down.F.sum + nDEGs_down.F
  nDEGs_down.M.sum <- nDEGs_down.M.sum + nDEGs_down.M
}
sex_DEGs.df.up <- data.frame(DEGs="Up-regulated",Cell_type="other",log2diffnDEGs=log2(nDEGs_up.F.sum/nDEGs_up.M.sum))
sex_DEGs.df.down <- data.frame(DEGs="Down-regulated",Cell_type="other",log2diffnDEGs=log2(nDEGs_down.F.sum/nDEGs_down.M.sum))

sex_DEGs.df <- rbind(sex_DEGs.df,sex_DEGs.df.up,sex_DEGs.df.down)


#devtools::install_github("kassambara/easyGgplot2")
library(easyGgplot2)

##      DEGs   Cell_type log2(Female DEGs/Male DEGs)
## 1 Up-regulated  Lunch      13.53
## 2 Female Dinner      16.81
## 3   Male  Lunch      16.24
## 4   Male Dinner      17.42

sex_DEGs.df$DEGs <- factor(sex_DEGs.df$DEGs,levels = c("Up-regulated","Down-regulated"))

ggplot2.barplot(data=sex_DEGs.df, xName='Cell_type', yName="log2diffnDEGs",
                groupName='DEGs', groupColors=c('#D55E00','#0072B2'),
                position=position_dodge(),
                #background and line colors
                backgroundColor="white", color="black", 
                xtitle="Cell Type", ytitle="log2(Female DEGs / Male DEGs)", 
                mainTitle="Sex difference in DEGs of stimulated PBMC",
                removePanelGrid=TRUE,removePanelBorder=TRUE,
                axisLine=c(0.5, "solid", "black")
) + theme(axis.text.x = element_text(size=9, angle=45, vjust = 1, hjust = 1))
ggsave("~/scRNAseq_analysis/nDEGs_F_M_stim.pdf", width = 7, height = 5,limitsize = FALSE)


### volcano plot ###
library(EnhancedVolcano)
markers.uGvs1G.stim.all.F.notaxid<-markers.uGvs1G.stim.all.F[grep("taxid|MT-",row.names(markers.uGvs1G.stim.all.F),invert = T),]
EnhancedVolcano(markers.uGvs1G.stim.all.F.notaxid,
                lab = rownames(markers.uGvs1G.stim.all.F.notaxid),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                ylab = '-log10(adj.p)',
                legendLabels = c("NS", expression(Log[2] ~ FC), "adj.p-value", expression(adj.p - value ~ and
                                                                                          ~ log[2] ~ FC)),
                xlim =c(-1.5,2),
                FCcutoff = 0.25,
                labSize = 5, 
                drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                title = "Female uG vs 1G",
                subtitle = bquote(italic("stimulated PBMC")))
ggsave("Female_uGvs1G_stimulated_Volcano.pdf",width=10,height=10)

markers.uGvs1G.stim.all.M.notaxid<-markers.uGvs1G.stim.all.M[grep("taxid|MT-",row.names(markers.uGvs1G.stim.all.M),invert = T),]
EnhancedVolcano(markers.uGvs1G.stim.all.M.notaxid,
                lab = rownames(markers.uGvs1G.stim.all.M.notaxid),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                ylab = '-log10(adj.p)',
                legendLabels = c("NS", expression(Log[2] ~ FC), "adj.p-value", expression(adj.p - value ~ and
                                                                                          ~ log[2] ~ FC)),
                xlim =c(-1.5,2),
                FCcutoff = 0.25,
                labSize = 5, 
                drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                title = "Male uG vs 1G",
                subtitle = bquote(italic("stimulated PBMC")))
ggsave("Male_uGvs1G_stimulated_Volcano.pdf",width=10,height=10)

