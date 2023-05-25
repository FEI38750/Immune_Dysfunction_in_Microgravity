library(Seurat)
library(patchwork)
library(sctransform)
library(glmGamPoi)
library(future)
library(dplyr)
library(ggplot2)

setwd("/opt/home/buckcenter.org/fwu/scRNAseq_analysis")

# change the current plan to access parallelization
plan("multisession", workers = 4)
options(future.globals.maxSize= +Inf) # increase maxSize of future function

Idents(Sample.combined)<-"orig.ident"
# get a table for all the genes in uG vs 1G unstimulated PBMCs in Female
markers.uGvs1G.unstim.all.F <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="F.uG_unstimulated", ident.2 ="F.1G_unstimulated",
                                           min.pct=0, logfc.threshold = 0,
                                           test.use = "MAST")
markers.uGvs1G.unstim.all.F<-markers.uGvs1G.unstim.all.F[order(markers.uGvs1G.unstim.all.F$avg_log2FC, decreasing = T),] # sort by log2FC

# save the marker lists of Female
# write.csv(markers.uGvs1G.unstim.F, "Female_markers_uGvs1G_unstimulated.csv")
write.csv(markers.uGvs1G.unstim.all.F, "Female_markers_uGvs1G_24K_unstimulated.csv")

# get a table for all the genes in uG vs 1G unstimulated PBMCs in Male
markers.uGvs1G.unstim.all.M <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="M.uG_unstimulated", ident.2 ="M.1G_unstimulated",
                                           min.pct=0, logfc.threshold = 0,
                                           test.use = "MAST")
markers.uGvs1G.unstim.all.M<-markers.uGvs1G.unstim.all.M[order(markers.uGvs1G.unstim.all.M$avg_log2FC, decreasing = T),] # sort by log2FC

# save the marker lists of Female
# write.csv(markers.uGvs1G.unstim.M, "Female_markers_uGvs1G_unstimulated.csv")
write.csv(markers.uGvs1G.unstim.all.M, "Male_markers_uGvs1G_24K_unstimulated.csv")


# DEGs for each cell types
Sample.combined$cell_sex <- paste0(Sample.combined$predicted.celltype.l2,"_",Sample.combined$orig.ident)
Idents(Sample.combined) <- "cell_sex"
cell4DEG<-unique(Sample.combined.unstim$predicted.celltype.l2)
for (d in cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating","pDC","Plasmablast")]){
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_F.uG_unstimulated"), ident.2 =paste0(d,"_F.1G_unstimulated"))
  write.csv(DEG.uG, paste0(d,"_uGvs1G_Female.csv"))
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_M.uG_unstimulated"), ident.2 =paste0(d,"_M.1G_unstimulated"))
  write.csv(DEG.uG, paste0(d,"_uGvs1G_Male.csv"))
}

# full gene list
for (d in cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating","pDC","Plasmablast")]){
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_F.uG_unstimulated"), ident.2 =paste0(d,"_F.1G_unstimulated"),
                        min.pct=0, logfc.threshold = 0,
                        test.use = "MAST")
  write.csv(DEG.uG, paste0(d,"_uGvs1G_Female.csv"))
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_M.uG_unstimulated"), ident.2 =paste0(d,"_M.1G_unstimulated"),
                        min.pct=0, logfc.threshold = 0,
                        test.use = "MAST")
  write.csv(DEG.uG, paste0(d,"_uGvs1G_Male.csv"))
}


# compare the Female DEGs/Male DEGs
setwd("../Sex_cellTypes_unstim") # can be change for stim
celltype4sex.unstim <- cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating","pDC","Plasmablast",
                                                 "MAIT","dnT","gdT","Treg","HSPC","NK_CD56bright","CD4 CTL","cDC2")]
celltype4sex.unstim <- append("Overall",celltype4sex.unstim)
sex_DEGs.df <- data.frame()
for (d in celltype4sex.unstim){
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

celltype4sex.unstim.other <- c("MAIT","dnT","gdT","Treg","HSPC","NK_CD56bright","CD4 CTL","cDC2")
nDEGs_up.F.sum <- 0
nDEGs_up.M.sum <- 0
nDEGs_down.F.sum <- 0
nDEGs_down.M.sum <- 0
for (o in celltype4sex.unstim.other){
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

sex_DEGs.df$DEGs <- factor(sex_DEGs.df$DEGs,levels = c("Up-regulated","Down-regulated"))

ggplot2.barplot(data=sex_DEGs.df, xName='Cell_type', yName="log2diffnDEGs",
                groupName='DEGs', groupColors=c('#D55E00','#0072B2'),
                position=position_dodge(),
                #background and line colors
                backgroundColor="white", color="black", 
                xtitle="Cell Type", ytitle="log2(Female DEGs / Male DEGs)", 
                mainTitle="Sex difference in DEGs of unstimulated PBMC",
                removePanelGrid=TRUE,removePanelBorder=TRUE,
                axisLine=c(0.5, "solid", "black")
) + theme(axis.text.x = element_text(size=9, angle=45, vjust = 1, hjust = 1))
ggsave("~/scRNAseq_analysis/nDEGs_F_M_unstim.pdf", width = 7.5, height = 5,limitsize = FALSE)




