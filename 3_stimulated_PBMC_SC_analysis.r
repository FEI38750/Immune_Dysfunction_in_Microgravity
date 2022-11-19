# ------ 3: stimulated PBMC SC analysis ------ #
# continue from the step 2_unstimulated_PBMC_SC_analysis.r
setwd("~/scRNAseq_analysis/stimulated")
# find DEGs/markers
## uG_stimulated vs 1G_stimulated
Idents(Sample.combined)<-"treatment"
markers.uG.stim <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="uG_stimulated", ident.2 ="1G_stimulated",
                               min.pct=0.005, logfc.threshold = 0.2,
                               test.use = "MAST")
markers.uG.stim<-markers.uG.stim[order(markers.uG.stim$avg_log2FC, decreasing = T),] # sort by log2FC
write.csv(markers.uG.stim, "markers_uGvs1G_stimulated.csv")

# Find all DEGs: get a table for all the genes in uG vs 1G stimulated PBMCs
## uG_stimulated vs 1G_stimulated
Idents(Sample.combined)<-"treatment"
markers.uGvs1G.stim.all <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="uG_stimulated", ident.2 ="1G_stimulated",
                               min.pct=0, logfc.threshold = 0,
                               test.use = "MAST")
markers.uGvs1G.stim.all<-markers.uGvs1G.stim.all[order(markers.uGvs1G.stim.all$avg_log2FC, decreasing = T),] # sort by log2FC
write.csv(markers.uGvs1G.stim.all, "markers_uGvs1G_stimulated_all.csv")

# cell proportion changes - stimulated
sub.prop2<-sub.prop.all[sub.prop.all$sample %in% c("1G_stimulated","uG_stimulated"),]
## prep cumulative sums for line links
cell.stack.wide2<-sub.prop2 %>%
  pivot_wider(names_from=sample, values_from=Freq) %>% 
  arrange(by=desc(Var1)) %>% 
  replace(is.na(.), 0) %>% # remove NA
  mutate(y=cumsum(`1G_stimulated`),
         yend=cumsum(uG_stimulated))
## set the levels in order we want
sub.prop2$sample<-factor(sub.prop2$sample, 
                         levels=c("1G_stimulated","uG_stimulated"))
sub.prop2$labs<-round(sub.prop2$Freq,3) # prepare cell proportion label
sub.prop2$labs[sub.prop2$labs<0.05]<-""
ggplot(sub.prop2, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(stat="identity",width=0.5, col="black") +
  geom_segment(data = cell.stack.wide2, aes(x=1.25, xend=1.75, y=y, yend=yend))+
  geom_text(aes(label= labs), size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="PBMC") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("stacked_barplot_stim.pdf",width=7,height=6)

# cell type ratio change: stimulated
ggplot(sub.prop2, aes(x = Var1, y = Freq, fill = sample)) +
  geom_bar(stat = "identity",position = "dodge") +
  ggtitle("Cell type frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(x="Cell types", y="Proportion", title="Cell type frequency") +
  guides(fill=guide_legend(title='Groups'))
ggsave("Cell_type_frequency_stimu.pdf",width=10,height=6)

# cell type ratio change: log2FC comparision stim; remove celltype ratio <1%
sub.prop2.1<-sub.prop2[sub.prop2$sample=="1G_stimulated",]
sub.prop2.1<-sub.prop2.1[sub.prop2.1$Freq>0.01,]
sub.prop2.2<-sub.prop2[sub.prop2$sample=="uG_stimulated",]
sub.prop2.2<-sub.prop2.2[sub.prop2.2$Freq>0.01,]
sub.prop.comm<-intersect(sub.prop2.1$Var1,sub.prop2.2$Var1)
sub.prop2.1<-sub.prop2.1[sub.prop2.1$Var1 %in% sub.prop.comm,]
sub.prop2.2<-sub.prop2.2[sub.prop2.2$Var1 %in% sub.prop.comm,]
sub.prop2vs1<-data.frame(Cell_type=sub.prop2.1$Var1,
                         group="uG vs 1G stimulated",
                         Fold_change=sub.prop2.2$Freq/sub.prop2.1$Freq)
sub.prop2vs1$log2FC<-log2(sub.prop2vs1$Fold_change)
ggplot(sub.prop2vs1, aes(Cell_type,log2FC))+
  geom_bar(position="dodge", stat="identity",fill="#999999") +
  ggtitle("uG vs 1G stimulated") +
  geom_text(aes(label=round(log2FC,2)),size=2.5,
            position=position_dodge(0.9), vjust=-0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust =0.5),
        axis.text.x = element_text(angle=45, hjust=1))
ggsave("CellType_ratio_uG_vs_1G_stimulated.pdf",width=6,height=4)

# subset stimulated samples
Idents(Sample.combined)<-"treatment"
Sample.combined.stim <- subset(Sample.combined, idents = c("1G_stimulated","uG_stimulated"))
# update the levels of treatment
Sample.combined.stim$treatment<-factor(Sample.combined.stim$treatment, 
                                       levels=c("1G_stimulated","uG_stimulated"))
DefaultAssay(Sample.combined.stim)
## view&save total umap in stimulated subset
Idents(Sample.combined.stim) <- "predicted.celltype.l2"
DimPlot(Sample.combined.stim,split.by = "treatment", label =T, repel = T)
ggsave("umap_predicted.celltype.stim.pdf",width=12,height=8)

#### draw dot plot for markers by Log2FC & p-value values
setwd("/opt/home/buckcenter.org/fwu/scRNAseq_analysis/stimulated")
features4plot.stim <- c("FTL","NCF1","CCL8","FTH1","CCL2","THBS1","MMP9","CTSB","IFI30","IDO1","TYROBP","CCL4","C15ORF48","SOD2","MARCKS",
                   "CTSZ","MS4A7","CCL7","DOCK4","CSTB","LGALS3","HLA-DRB1","CXCL8","ANXA5","PLA2G7","IL1B","PSAP","CCL4L2","S100A9",
                   "CXCL10","PABPC1",
                   "CLEC2B","MT-CO3","HLA-C","RPS4Y1","LAG3","HAPLN3","CD38","TAP1","SAT1","GSTK1","AL138963.4","APOL6","MTRNR2L8",
                   "RGS1","TNFSF10","SLFN5","HBB","RBM3","WARS","XRN1","SYNE2","SOCS1","IRF1","PLAAT4","PRF1","GBP2","GBP1","GBP4",
                   "GBP5")
celltype4dot.stim<-c("CD4 TCM", "CD4 Naive", "CD4 TEM", "B naive","B intermediate",
                       "B memory","CD8 Naive","CD8 TEM", "CD8 TCM", "dnT", "MAIT",
                       "cDC2", "pDC", "NK", "NK_CD56bright", "HSPC", "Plasmablast", "CD14 Mono", "CD16 Mono")
celltype4dot.stim<-append("Overall",celltype4dot.stim)
lsmarker.dot.stim<-data.frame()
for (d in celltype4dot.stim){
  markers4dot.stim <- read.csv(paste0(d,"_uGvs1G_stim.csv"))
  if (nrow(markers4dot.stim) != 0){
    genes<-filter(markers4dot.stim,X %in% features4plot.stim)
    genes<-genes[order(genes$avg_log2FC,decreasing = T),]
    genes$Cell_type<-d
    lsmarker.dot.stim<-rbind(lsmarker.dot.stim,genes)
  }
}

# replace the 0 p-values with lowest number
lsmarker.dot.stim$p_val_adj[lsmarker.dot.stim$p_val_adj==0] <- sort(unique(lsmarker.dot.stim$p_val_adj))[2]

gs.stim<-lsmarker.dot.stim$X %>% unique()
lsmarker.dot.stim$Cell_type<-factor(lsmarker.dot.stim$Cell_type, levels = unique(lsmarker.dot.stim$Cell_type)) # order the x axis

lsmarker.dot.stim %>% filter(X %in% gs.stim) %>% filter(p_val_adj <0.05) %>%
  ggplot(aes(x=X, y = Cell_type, color = avg_log2FC, size = -log10(p_val_adj))) + 
  geom_point() +
  scale_y_discrete(labels=function(y) str_wrap(y, width=70)) +
  ylab('Cell Types') +
  xlab('Genes') +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size=9, angle=45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=10))+
  theme(axis.line  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  scale_size_area(
    max_size = 5,
    breaks = c(0,25,50,100),
    labels = c("0","25","50","100+"),
    guide = "legend",
    limits = c(0, 100),
    oob = scales::squish
  )+
  colorspace::scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-1,1), oob = scales::squish, name = 'Log2FC',
                                                na.value="transparent")
ggsave("Dot_Markers_uGvs1G_stim.pdf", width = 14, height = 5,limitsize = FALSE)

# ------ find DEG for all cell types during gravity change ------ #
# with default filter min.pct=0.1 and logfc=0.25
Idents(Sample.combined)<-"celltype.treatment"
plan("multisession", workers = 20)
table(Sample.combined.stim$predicted.celltype.l2) # check cell numbers in each types
cell4DEG.stim<-unique(Sample.combined.stim$predicted.celltype.l2)
for (d in cell4DEG.stim[!cell4DEG.stim %in% c("Unknow","cDC1","CD4 Proliferating","NK Proliferating")]){
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_uG_stimulated"), ident.2 =paste0(d,"_1G_stimulated"))
  write.csv(DEG.uG, paste0(d,"_uGvs1G_stim.csv"))
}
# include all the genes
setwd("~/scRNAseq_analysis/Markers_cellType_stim_all")
for (d in cell4DEG.stim[!cell4DEG.stim %in% c("Unknow","cDC1","CD4 Proliferating","NK Proliferating")]){
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_uG_stimulated"), ident.2 =paste0(d,"_1G_stimulated"),
                        min.pct = 0,logfc.threshold =0)
  write.csv(DEG.uG, paste0(d,"_uGvs1G_stim_allgenes.csv"))
}

# volcano plot for overall uG vs 1G stim
library(EnhancedVolcano)
markers.uGvs1G.stim.all.notaxid<-markers.uGvs1G.unstim.all[grep("taxid|MT-",row.names(markers.uGvs1G.stim.all),invert = T),]

EnhancedVolcano(markers.uGvs1G.stim.all.notaxid,
                lab = rownames(markers.uGvs1G.stim.all.notaxid),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim =c(-1.5,2),
                FCcutoff = 0.25,
                labSize = 5, 
                drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                title = "uG vs 1G",
                subtitle = bquote(italic("unstimulated PBMC")))

ggsave("uGvs1G_stimulated_Volcano.pdf",width=10,height=10)


### --- iAge index calculation --- ###
library(scales)
library(ggpubr)
# read GE_iAge gene and coefficients
GE_iAge <- read.csv("GE_iAge.csv")

# diveide into to groups: 1G and uG
stim_1G4iAge<-Sample.combined.stim@assays[["RNA"]]@data[intersect(GE_iAge$Genes,row.names(Sample.combined.stim@assays[["RNA"]]@data)),Sample.combined.stim$gravity=="1G"]
stim_uG4iAge<-Sample.combined.stim@assays[["RNA"]]@data[intersect(GE_iAge$Genes,row.names(Sample.combined.stim@assays[["RNA"]]@data)),Sample.combined.stim$gravity=="uG"]

# calculate iAge index for each cell
sumGene<-0
for (g in intersect(GE_iAge$Genes,row.names(Sample.combined.stim@assays[["RNA"]]@data))){
  resGene<-rescale(stim_1G4iAge[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
  sumGene<-sumGene+resGene
}
iAge_index_1G<-sumGene+10
summary(iAge_index_1G)

sumGene<-0
for (g in intersect(GE_iAge$Genes,row.names(Sample.combined.stim@assays[["RNA"]]@data))){
  resGene<-rescale(stim_uG4iAge[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
  sumGene<-sumGene+resGene
}
iAge_index_uG<-sumGene+10
summary(iAge_index_uG)

iAge.data<-rbind(data.frame(iAge_index=iAge_index_1G,group=rep("1G",length(iAge_index_1G))),
                 data.frame(iAge_index=iAge_index_uG,group=rep("uG",length(iAge_index_uG))))
# plot&save iAge index for overall/general cells in 1G and uG
ggboxplot(iAge.data, x = "group", y = "iAge_index",
          color = "group", palette = "jco",outlier.shape = NA) +
  stat_compare_means (method = "wilcox.test", label = "p.signif")
ggsave("GE_iAge_PBMC_stim.pdf", width = 6, height = 5)
write.csv(iAge.data,"GE_iAge_PBMC_stim.csv")

# --- iAge for cell types --- #
# calculate iAge index for each cells
iAge.data.cell<- data.frame(matrix(ncol = 3, nrow = 0))
colnames(iAge.data.cell) <- c("iAge_index","cell_type","group")
for (c in cell4gsea.stim[-1]){
  cell4iage<-dplyr::filter(Sample.combined.stim@meta.data,predicted.celltype.l2==c & gravity == "1G")
  cell4iage<-row.names(cell4iage)
  stim_1G4iAge_cell<-stim_1G4iAge[,cell4iage]
  sumGene<-0
  for (g in intersect(GE_iAge$Genes,row.names(Sample.combined.stim@assays[["RNA"]]@data))){
    resGene<-rescale(stim_1G4iAge_cell[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
    sumGene<-sumGene+resGene
  }
  iAge_index_1G<-sumGene+10
  iAge.data.cell<-rbind(iAge.data.cell,
                        data.frame(iAge_index=iAge_index_1G,cell_type=rep(c,length(iAge_index_1G)),group=rep("1G",length(iAge_index_1G))))
}

iAge.data.cell2<- data.frame(matrix(ncol = 3, nrow = 0))
colnames(iAge.data.cell2) <- c("iAge_index","cell_type","group")
for (c in cell4gsea.stim[-1]){
  cell4iage<-dplyr::filter(Sample.combined.stim@meta.data,predicted.celltype.l2==c & gravity == "uG")
  cell4iage<-row.names(cell4iage)
  stim_uG4iAge_cell<-stim_uG4iAge[,cell4iage]
  sumGene<-0
  for (g in intersect(GE_iAge$Genes,row.names(Sample.combined.stim@assays[["RNA"]]@data))){
    resGene<-rescale(stim_uG4iAge_cell[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
    sumGene<-sumGene+resGene
  }
  iAge_index_uG<-sumGene+10
  iAge.data.cell2<-rbind(iAge.data.cell2,
                         data.frame(iAge_index=iAge_index_uG,cell_type=rep(c,length(iAge_index_uG)),group=rep("uG",length(iAge_index_uG))))
}
iAge.data.cell<-rbind(iAge.data.cell,iAge.data.cell2)
write.csv(iAge.data.cell,"GE_iAge_PBMC_celltypes_stim.csv")
# plot&save iAge index for each cell type in 1G and uG
ggboxplot(iAge.data.cell, x = "cell_type", y = "iAge_index",
          color = "group", palette = "jco",outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("GEiage_PBMC_celltypes_stim.pdf", width = 8, height = 6)


# ------ Monocle3 single-cell trajectory analysis ------ #
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

# uG vs 1G in Stim
Idents(Sample.combined.stim) <- "gravity"
DefaultAssay(Sample.combined.stim)<-"integrated"

cds <- as.cell_data_set(subset(Sample.combined.stim, idents = "uG"))
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)
ggsave("PBMC_trej_stim_uG.pdf", width =7, height = 6)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE)
ggsave("PBMC_trej_stim_uG_sudotime.pdf", width =7, height = 6)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
ggsave("PBMC_trej_stim_uG_partition.pdf", width =7, height = 6)

cds <- as.cell_data_set(subset(Sample.combined.stim, idents = "1G"))
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)
ggsave("PBMC_trej_stim_1G.pdf", width =7, height = 6)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE)
ggsave("PBMC_trej_stim_1G_sudotime.pdf", width =7, height = 6)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
ggsave("PBMC_trej_stim_1G_partition.pdf", width =7, height = 6)

### --- Cellular senescence score calculation by using SenMayo gene set --- ###
### stim
pbmc_module_stim <- AddModuleScore(
  object = Sample.combined.stim,
  features = geneset_features,
  name = 'SenMayo')
ggplot(pbmc_module_stim@meta.data, aes(x=SenMayo1, fill=gravity)) +
  geom_density(alpha=0.4) + theme_bw()
ggsave("PBMC_stim_SenMayo.pdf", width = 6, height = 5)
ggplot(pbmc_module_stim@meta.data) + 
  geom_density_ridges(aes(x=SenMayo1, y= predicted.celltype.l2, fill=gravity), alpha=0.5) + theme_bw()
ggsave("PBMC_stim_SenMayo_celltypes.pdf", width = 6, height = 6)
t.test(filter(pbmc_module_stim@meta.data,gravity=="1G")$SenMayo1,
       filter(pbmc_module_stim@meta.data,gravity=="uG")$SenMayo1)

c1<-c(rep("1G",length(filter(pbmc_module_stim@meta.data,gravity=="1G")$SenMayo1)),
      rep("uG",length(filter(pbmc_module_stim@meta.data,gravity=="uG")$SenMayo1)))
c2<-c(filter(pbmc_module_stim@meta.data,gravity=="1G")$SenMayo1,filter(pbmc_module_stim@meta.data,gravity=="uG")$SenMayo1)
dat<-data.frame(Group=c1,Score=c2)
ggplot(dat,aes(Group,Score)) + 
  geom_violin(aes(fill = Group)) +
  theme_bw() + scale_fill_brewer(palette = 'Set2',direction = -1) +
  labs(title="SenMayo") +
  geom_boxplot(width=0.1,outlier.shape = NA) +
  ggpubr::stat_compare_means(method = "t.test", 
                             comparisons=list(c("1G","uG")), 
                             label = 'p.signif')
ggsave("uGvs1G_stim_SenMayo_violin.pdf", width = 5, height = 7)

ggboxplot(dat, x="Group",y="Score",
          color = "Group", palette = "jco",outlier.shape = NA) + coord_cartesian(ylim = c(-0.1,0.3)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",label.y=0.25,
                     comparisons=list(c("1G","uG")))
ggsave("uGvs1G_stim_SenMayo.pdf", width = 6, height = 5)