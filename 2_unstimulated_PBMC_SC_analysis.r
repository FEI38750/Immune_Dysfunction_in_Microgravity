# ------ 2: unstimulated PBMC SC analysis ------ #
# continue from the step 1_Preprocess_the_count_matrix.r
setwd("~/scRNAseq_analysis/unstimulated")
# find DEGs/markers
Idents(Sample.combined)<-"predicted.celltype.l2"
# Prepare object to run differential expression on SCT assay with multiple models
Sample.combined <- PrepSCTFindMarkers(Sample.combined)
# find markers for every cluster compared to all remaining cells
plan("multisession", workers = 16)
options(future.globals.maxSize= +Inf) # increase maxSize of future function
# find markers for each predicted cell type
markers.all.celltype <- FindAllMarkers(Sample.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(markers.all.celltype, "~/scRNAseq_analysis/markers_all_celltype.csv")

Idents(Sample.combined)<-"treatment"
# gene DEGs/markers of uG vs 1G unstimulated PBMCs
markers.uGvs1G.unstim <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="uG_unstimulated", ident.2 ="1G_unstimulated",
                               min.pct=0.005, logfc.threshold = 0.1,
                               test.use = "MAST")
markers.uGvs1G.unstim<-markers.uGvs1G.unstim[order(markers.uGvs1G.unstim$avg_log2FC, decreasing = T),] # sort by log2FC

# get a table for all the genes in uG vs 1G unstimulated PBMCs
markers.uGvs1G.unstim.all <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="uG_unstimulated", ident.2 ="1G_unstimulated",
                               min.pct=0, logfc.threshold = 0,
                               test.use = "MAST")
markers.uGvs1G.unstim.all<-markers.uGvs1G.unstim.all[order(markers.uGvs1G.unstim.all$avg_log2FC, decreasing = T),] # sort by log2FC

# save the marker lists
write.csv(markers.uGvs1G.unstim, "markers_uGvs1G_390_unstimulated.csv")
write.csv(markers.uGvs1G.unstim.all, "markers_uGvs1G_24K_unstimulated.csv")

# save gene counts for iAge calculation
write.csv(as.data.frame(Sample.combined@assays[["SCT"]]@data),"PBMC_SCT_data.csv")
write.csv(as.data.frame(Sample.combined@meta.data),"PBMC_metadata.csv")

# ------ Cell type ratio change stack barplot ------ #
sub.prop.all<-data.frame()
for (l in levels(Sample.combined$treatment)){
  sub.treat<-Sample.combined@meta.data[Sample.combined$treatment==l,]
  sub.prop<-data.frame(table(sub.treat$predicted.celltype.l2)/sum(table(sub.treat$predicted.celltype.l2)))
  sub.prop$sample<-l
  sub.prop.all<-rbind(sub.prop.all,sub.prop)
}

library(tidyr)
# All sample groups: cell type ratio change stack barplot
## set the levels in order we want
sub.prop.all$sample<-factor(sub.prop.all$sample, 
                                  levels=c("1G_unstimulated","uG_unstimulated","1G_stimulated","uG_stimulated"))
ggplot(sub.prop.all, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="PBMC") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("stacked_barplot_all.pdf",width=8,height=6)

# cell proportion changes - unstimulated
sub.prop1<-sub.prop.all[sub.prop.all$sample %in% c("1G_unstimulated","uG_unstimulated"),]
## prep cumulative sums for line links
cell.stack.wide1<-sub.prop1 %>%
  pivot_wider(names_from=sample, values_from=Freq) %>% 
  arrange(by=desc(Var1)) %>% 
  replace(is.na(.), 0) %>% # remove NA
  mutate(y=cumsum(`1G_unstimulated`),
         yend=cumsum(uG_unstimulated))
## set the levels in order we want
sub.prop1$sample<-factor(sub.prop1$sample, 
                            levels=c("1G_unstimulated","uG_unstimulated"))
sub.prop1$labs<-round(sub.prop1$Freq,3) # prepare cell proportion label
sub.prop1$labs[sub.prop1$labs<0.05]<-""
ggplot(sub.prop1, aes(fill=Var1, y=Freq, x=sample)) + 
  geom_bar(stat="identity",width=0.5, col="black") +
  geom_segment(data = cell.stack.wide1, aes(x=1.25, xend=1.75, y=y, yend=yend))+
  geom_text(aes(label= labs), size = 3, position = position_stack(vjust = 0.5)) +
  theme_bw() +
  labs(x="Sample", y="Proportion", title="PBMC") +
  guides(fill=guide_legend(title='Cell types'))
ggsave("stacked_barplot_unstim.pdf",width=7,height=6)

# cell type ratio change: unstimulated
ggplot(sub.prop1, aes(x = Var1, y = Freq, fill = sample)) +
  geom_bar(stat = "identity",position = "dodge") +
  ggtitle("Cell type frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(x="Cell types", y="Proportion", title="Cell type frequency") +
  guides(fill=guide_legend(title='Groups'))
ggsave("Cell_type_frequency_unstimu.pdf",width=10,height=6)

# cell type ratio change: log2FC comparision unstim; remove celltype ratio <1%
sub.prop1.1<-sub.prop1[sub.prop1$sample=="1G_unstimulated",]
sub.prop1.1<-sub.prop1.1[sub.prop1.1$Freq>0.01,]
sub.prop1.2<-sub.prop1[sub.prop1$sample=="uG_unstimulated",]
sub.prop1.2<-sub.prop1.2[sub.prop1.2$Freq>0.01,]
sub.prop.comm<-intersect(sub.prop1.1$Var1,sub.prop1.2$Var1)
sub.prop1.1<-sub.prop1.1[sub.prop1.1$Var1 %in% sub.prop.comm,]
sub.prop1.2<-sub.prop1.2[sub.prop1.2$Var1 %in% sub.prop.comm,]
sub.prop2vs1<-data.frame(Cell_type=sub.prop1.1$Var1,
                         group="uG vs 1G unstimulated",
                         Fold_change=sub.prop1.2$Freq/sub.prop1.1$Freq)
sub.prop2vs1$log2FC<-log2(sub.prop2vs1$Fold_change)
ggplot(sub.prop2vs1, aes(Cell_type,log2FC))+
  geom_bar(position="dodge", stat="identity",fill="#999999") +
  ggtitle("uG vs 1G unstimulated") +
  geom_text(aes(label=round(log2FC,2)),size=2.5,
            position=position_dodge(0.9), vjust=-0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust =0.5),
        axis.text.x = element_text(angle=45, hjust=1))
ggsave("CellType_ratio_uG_vs_1G_unstimulated.pdf",width=6,height=4)

# subset unstimulated samples
Idents(Sample.combined)<-"treatment"
Sample.combined.unstim <- subset(Sample.combined, idents = c("1G_unstimulated","uG_unstimulated"))
# update the levels of treatment
Sample.combined.unstim$treatment<-factor(Sample.combined.unstim$treatment, 
                            levels=c("1G_unstimulated","uG_unstimulated"))
## view&save total umap in unstimulated subset
Idents(Sample.combined.unstim) <- "predicted.celltype.l2"
DimPlot(Sample.combined.unstim,split.by = "treatment", label =T, repel = T)
ggsave("umap_predicted.celltype.unstim.pdf",width=12,height=8)

#### Dot plot:
##### create a new column of meta data (1G, uG) for aesthetic and more...
Sample.combined.unstim$gravity<-with(Sample.combined.unstim@meta.data, ifelse(treatment=="uG_unstimulated","uG","1G"))
#### draw dot plot for markers by Log2FC & p-value values
features4plot <- c("MMP9","THBS1","CCL2","FTH1","S100A9","FTL","S100A8","SOD2","MS4A7","C15ORF48","CYP27A1",
             "CYP1B1","CYBB","TMEM176B","PLA2G7","SLC7A11","NFKBIA","FCN1","TALDO1","ANXA5","S100A12",
             "TXN","CALR","SERF2","HSP90AB1","TMEM176A","MARCKS","KYNU","MRC1","IFI6","RPL38","FCER1G",
             "CCL4","HSP90AA1","IER3","TGFBI","HINT1","ATP5F1E","PSME2","PABPC1","GBP1","CIRBP","RPL3",
             "HNRNPA1","RBM3","MTRNR2L8","STAT1")
celltype4dot.unstim<-c("CD4 TCM", "CD4 Naive", "CD4 TEM", "CD4 CTL", "B naive","B intermediate",
                   "B memory","CD8 Naive","CD8 TEM", "CD8 TCM", "Treg", "gdT", "dnT", "MAIT",
                   "cDC2", "pDC", "NK", "NK_CD56bright", "HSPC", "Plasmablast", "CD14 Mono", "CD16 Mono")
celltype4dot.unstim<-append("Overall",celltype4dot.unstim)
lsmarker.dot<-data.frame()
for (d in celltype4dot.unstim){
  markers4dot <- read.csv(paste0(d,"_uGvs1G.csv"))
  if (nrow(markers4dot) != 0){
    genes<-filter(markers4dot,X %in% features4plot)
    genes<-genes[order(genes$avg_log2FC,decreasing = T),]
    genes$Cell_type<-d
    lsmarker.dot<-rbind(lsmarker.dot,genes)
  }
}
# replace the 0 p-values with lowest number
lsmarker.dot$p_val_adj[lsmarker.dot$p_val_adj==0] <- sort(unique(lsmarker.dot$p_val_adj))[2]

gs<-lsmarker.dot$X %>% unique()
lsmarker.dot$Cell_type<-factor(lsmarker.dot$Cell_type, levels = unique(lsmarker.dot$Cell_type)) # order the x axis

lsmarker.dot %>% filter(X %in% gs) %>% filter(p_val_adj <0.05) %>%
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
ggsave("Dot_Markers_uGvs1G_unstim.pdf", width = 11, height = 6,limitsize = FALSE)
scale_color_viridis()

# ------ find DEG for all cell types during gravity change ------ #
# with default filter min.pct=0.1 and logfc=0.25
Idents(Sample.combined)<-"celltype.treatment"
plan("multisession", workers = 20)
cell4DEG<-unique(Sample.combined.unstim$predicted.celltype.l2)
for (d in cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating")]){
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_uG_unstimulated"), ident.2 =paste0(d,"_1G_unstimulated"))
  write.csv(DEG.uG, paste0(d,"_uGvs1G.csv"))
}
# include all the genes
setwd("~/scRNAseq_analysis/Markers_cellType_unstim_all")
for (d in cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating")]){
  DEG.uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_uG_unstimulated"), ident.2 =paste0(d,"_1G_unstimulated"),
                        min.pct = 0,logfc.threshold =0)
  write.csv(DEG.uG, paste0(d,"_uGvs1G_allgenes.csv"))
}

# volcano plot for overall uG vs 1G unstim
library(EnhancedVolcano)
markers.uGvs1G.unstim.all.notaxid<-markers.uGvs1G.unstim.all[grep("taxid|MT-",row.names(markers.uGvs1G.unstim.all),invert = T),]

EnhancedVolcano(markers.uGvs1G.unstim.all.notaxid,
                lab = rownames(markers.uGvs1G.unstim.all.notaxid),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlim =c(-1.5,2),
                FCcutoff = 0.25,
                labSize = 5, 
                drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                title = "uG vs 1G",
                subtitle = bquote(italic("unstimulated PBMC")))

ggsave("uGvs1G_unstimulated_Volcano.pdf",width=10,height=10)


## microbial abundance comparison
Mycob.1G <-Sample.combined.unstim[["SCT"]]@data["Mycobacterium canettii (taxid 78331)",row.names(Sample.combined.unstim@meta.data)[Sample.combined.unstim$gravity=="1G"]]
Mycob.uG <-Sample.combined.unstim[["SCT"]]@data["Mycobacterium canettii (taxid 78331)",row.names(Sample.combined.unstim@meta.data)[Sample.combined.unstim$gravity=="uG"]]
t.test(Mycob.1G,Mycob.uG) # quick check here; data can also be exported to Prism
c1<-c(rep(deparse(substitute(Mycob.1G)),length(Mycob.1G)),
      rep(deparse(substitute(Mycob.uG)),length(Mycob.uG)))
c2<-c(Mycob.1G,Mycob.uG)
dat<-data.frame(Group=c1,Express_level=c2)
ggplot(dat,aes(Group,Express_level)) + 
  geom_violin(aes(fill = Group)) +
  theme_bw() + scale_fill_brewer(palette = 'Set2') +
  labs(title="Mycobacterium canettii") +
  ggpubr::stat_compare_means(method = "t.test", 
                             comparisons=list(c(deparse(substitute(Mycob.1G)),deparse(substitute(Mycob.uG)))), 
                             label = 'p.signif', label.y=0.00095)
ggsave("Mycobacterium_canettii_comparison.pdf", width = 6, height = 6)

# whole microbiome; preliminary test
microbiome.unstim <- subset(as.data.frame(Sample.combined.unstim[["SCT"]]@data),row.names(Sample.combined.unstim[["SCT"]]@data) %in% grep("taxid", row.names(Sample.combined.unstim[["SCT"]]@data),value = T))
t.test(colSums(microbiome.unstim[names(cell.1G)]),
       colSums(microbiome.unstim[names(cell.uG)]))
# recognized all "virus"; preliminary test
microbiome.unstim.vir <- microbiome.unstim[grepl("virus",row.names(microbiome.unstim)),]
t.test(colSums(microbiome.unstim.vir[names(cell.1G)]),
       colSums(microbiome.unstim.vir[names(cell.uG)]))

# quick look abundance
head(sort(rowSums(microbiome.unstim),decreasing = T), 100) # all
head(sort(rowSums(microbiome.unstim.vir),decreasing = T), 20) # Top20 virus only

# "Gammaretrovirus (taxid 153135)": the most abundance virus in the samples; preliminary test
microbiome.unstim.gamma <- microbiome.unstim["Gammaretrovirus (taxid 153135)",]
t.test(colSums(microbiome.unstim.gamma[names(cell.1G)]),
       colSums(microbiome.unstim.gamma[names(cell.uG)]))

## --- microbiome abundance in fraction --- ##
microbiome.unstim.counts <- subset(as.data.frame(Sample.combined.unstim[["RNA"]]@counts),row.names(Sample.combined.unstim[["RNA"]]@counts) %in% grep("taxid", row.names(Sample.combined.unstim[["RNA"]]@counts),value = T))
microbiome.unstim.frac <- microbiome.unstim.counts/colSums(as.data.frame(Sample.combined.unstim[["RNA"]]@counts))

### "Gammaretrovirus (taxid 153135)"; comparison test; data can also be exported to Prism
microbiome.unstim.gamma.frac <- microbiome.unstim.frac["Gammaretrovirus (taxid 153135)",]
gamma.1G<-colSums(microbiome.unstim.gamma.frac[names(cell.1G)])
gamma.uG<-colSums(microbiome.unstim.gamma.frac[names(cell.uG)])
t.test(gamma.1G,gamma.uG)
write.csv(gamma.1G,"gamma.1G.csv"); write.csv(gamma.uG,"gamma.uG.csv")

### "Mycobacterium canettii (taxid 78331)"; comparison test; data can also be exported to Prism
microbiome.unstim.mycob.frac <- microbiome.unstim.frac["Mycobacterium canettii (taxid 78331)",]
mycob.1G<-colSums(microbiome.unstim.mycob.frac[names(cell.1G)])
mycob.uG<-colSums(microbiome.unstim.mycob.frac[names(cell.uG)])
t.test(mycob.1G,mycob.uG)
write.csv(mycob.1G,"mycob.1G.csv"); write.csv(mycob.uG,"mycob.uG.csv")

# recognized all "virus"; comparison test; data can also be exported to Prism
microbiome.unstim.vir.frac <- microbiome.unstim.frac[grepl("virus",row.names(microbiome.unstim.frac)),]
t.test(colSums(microbiome.unstim.vir.frac[names(cell.1G)]),
       colSums(microbiome.unstim.vir.frac[names(cell.uG)]))
gene.1G<-colSums(microbiome.unstim.vir.frac[names(cell.1G)])
gene.uG<-colSums(microbiome.unstim.vir.frac[names(cell.uG)])
t.test(gene.1G,gene.uG)
write.csv(gene.1G,"virus.1G.csv"); write.csv(gene.uG,"virus.uG.csv")


### --- iAge index calculation --- ###
library(scales)
library(ggpubr)
# read GE_iAge gene and coefficients
GE_iAge <- read.csv("GE_iAge.csv")

# diveide into to groups: 1G and uG
unstim_1G4iAge<-Sample.combined.unstim@assays[["RNA"]]@data[intersect(GE_iAge$Genes,row.names(Sample.combined.unstim@assays[["RNA"]]@data)),Sample.combined.ununstim$gravity=="1G"]
unstim_uG4iAge<-Sample.combined.unstim@assays[["RNA"]]@data[intersect(GE_iAge$Genes,row.names(Sample.combined.unstim@assays[["RNA"]]@data)),Sample.combined.unstim$gravity=="uG"]

# calculate iAge index for each cell
sumGene<-0
for (g in intersect(GE_iAge$Genes,row.names(Sample.combined.unstim@assays[["RNA"]]@data))){
  resGene<-rescale(unstim_1G4iAge[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
  sumGene<-sumGene+resGene
}
iAge_index_1G<-sumGene+10
summary(iAge_index_1G)

sumGene<-0
for (g in intersect(GE_iAge$Genes,row.names(Sample.combined.unstim@assays[["RNA"]]@data))){
  resGene<-rescale(unstim_uG4iAge[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
  sumGene<-sumGene+resGene
}
iAge_index_uG<-sumGene+10
summary(iAge_index_uG)

iAge.data <- rbind(data.frame(iAge_index=iAge_index_1G,group=rep("1G",length(iAge_index_1G))),
                 data.frame(iAge_index=iAge_index_uG,group=rep("uG",length(iAge_index_uG))))
# plot&save iAge index for overall/general cells in 1G and uG
ggboxplot(iAge.data, x = "group", y = "iAge_index",
          color = "group", palette = "jco",outlier.shape = NA) +
  stat_compare_means (method = "wilcox.test", label = "p.signif")
ggsave("GE_iAge_PBMC_unstim.pdf", width = 6, height = 5)
write.csv(iAge.data,"GE_iAge_PBMC_unstim.csv")

# --- iAge for cell types --- #
# calculate iAge index for each cells
iAge.data.cell<- data.frame(matrix(ncol = 3, nrow = 0))
colnames(iAge.data.cell) <- c("iAge_index","cell_type","group")
for (c in cell4gsea.unstim[-1]){
  cell4iage<-dplyr::filter(Sample.combined.unstim@meta.data,predicted.celltype.l2==c & gravity == "1G")
  cell4iage<-row.names(cell4iage)
  unstim_1G4iAge_cell<-unstim_1G4iAge[,cell4iage]
  sumGene<-0
  for (g in intersect(GE_iAge$Genes,row.names(Sample.combined.unstim@assays[["RNA"]]@data))){
    resGene<-rescale(unstim_1G4iAge_cell[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
    sumGene<-sumGene+resGene
  }
  iAge_index_1G<-sumGene+10
  iAge.data.cell<-rbind(iAge.data.cell,
                        data.frame(iAge_index=iAge_index_1G,cell_type=rep(c,length(iAge_index_1G)),group=rep("1G",length(iAge_index_1G))))
}

iAge.data.cell2<- data.frame(matrix(ncol = 3, nrow = 0))
colnames(iAge.data.cell2) <- c("iAge_index","cell_type","group")
for (c in cell4gsea.unstim[-1]){
  cell4iage<-dplyr::filter(Sample.combined.unstim@meta.data,predicted.celltype.l2==c & gravity == "uG")
  cell4iage<-row.names(cell4iage)
  unstim_uG4iAge_cell<-unstim_uG4iAge[,cell4iage]
  sumGene<-0
  for (g in intersect(GE_iAge$Genes,row.names(Sample.combined.unstim@assays[["RNA"]]@data))){
    resGene<-rescale(unstim_uG4iAge_cell[g,], to = c(0.01, 5)*GE_iAge[GE_iAge$Genes==g,"Coefficients"])
    sumGene<-sumGene+resGene
  }
  iAge_index_uG<-sumGene+10
  iAge.data.cell2<-rbind(iAge.data.cell2,
                         data.frame(iAge_index=iAge_index_uG,cell_type=rep(c,length(iAge_index_uG)),group=rep("uG",length(iAge_index_uG))))
}
iAge.data.cell<-rbind(iAge.data.cell,iAge.data.cell2)
write.csv(iAge.data.cell,"GE_iAge_PBMC_celltypes_unstim.csv")
# plot&save iAge index for each cell type in 1G and uG
ggboxplot(iAge.data.cell, x = "cell_type", y = "iAge_index",
          color = "group", palette = "jco",outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave("GEiage_PBMC_celltypes_unstim.pdf", width = 8, height = 6)

# ------ Monocle3 single-cell trajectory analysis ------ #
library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

Idents(Sample.combined.unstim) <- "gravity"
DefaultAssay(Sample.combined.unstim)<-"integrated"
# uG vs 1G in Unstim
cds <- as.cell_data_set(subset(Sample.combined.unstim, idents = "uG"))
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)
ggsave("PBMC_trej_unstim_uG.pdf", width =7, height = 6)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE)
ggsave("PBMC_trej_unstim_uG_sudotime.pdf", width =7, height = 6)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
ggsave("PBMC_trej_unstim_uG_partition.pdf", width =7, height = 6)

cds <- as.cell_data_set(subset(Sample.combined.unstim, idents = "1G"))
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds)
ggsave("PBMC_trej_unstim_1G.pdf", width =7, height = 6)
plot_cells(cds,color_cells_by = "pseudotime",label_cell_groups=FALSE)
ggsave("PBMC_trej_unstim_1G_sudotime.pdf", width =7, height = 6)
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
ggsave("PBMC_trej_unstim_1G_partition.pdf", width =7, height = 6)

### --- Cellular senescence score calculation by using SenMayo gene set --- ###
geneset_features <- list(c("ACVR1B","ANG","ANGPT1","ANGPTL4","AREG","AXL","BEX3","BMP2",
                           "BMP6","C3","CCL1","CCL13","CCL16","CCL2","CCL20","CCL24","CCL26",
                           "CCL3","CCL3L1","CCL4","CCL5","CCL7","CCL8","CD55","CD9","CSF1",
                           "CSF2","CSF2RB","CST4","CTNNB1","CTSB","CXCL1","CXCL10","CXCL12",
                           "CXCL16","CXCL2","CXCL3","CXCL8","CXCR2","DKK1","EDN1","EGF","EGFR",
                           "EREG","ESM1","ETS2","FAS","FGF1","FGF2","FGF7","GDF15","GEM","GMFG",
                           "HGF","HMGB1","ICAM1","ICAM3","IGF1","IGFBP1","IGFBP2","IGFBP3",
                           "IGFBP4","IGFBP5","IGFBP6","IGFBP7","IL10","IL13","IL15","IL18",
                           "IL1A","IL1B","IL2","IL32","IL6","IL6ST","IL7","INHA","IQGAP2","ITGA2",
                           "ITPKA","JUN","KITLG","LCP1","MIF","MMP1","MMP10","MMP12","MMP13",
                           "MMP14","MMP2","MMP3","MMP9","NAP1L4","NRG1","PAPPA","PECAM1","PGF",
                           "PIGF","PLAT","PLAU","PLAUR","PTBP1","PTGER2","PTGES","RPS6KA5","SCAMP4",
                           "SELPLG","SEMA3F","SERPINB4","SERPINE1","SERPINE2","SPP1","SPX","TIMP2",
                           "TNF","TNFRSF10C","TNFRSF11B","TNFRSF1A","TNFRSF1B","TUBGCP2","VEGFA",
                           "VEGFC","VGF","WNT16","WNT2"))
### unstim
pbmc_module <- AddModuleScore(
  object = Sample.combined.unstim,
  features = geneset_features,
  name = 'SenMayo')
ggplot(pbmc_module@meta.data, aes(x=SenMayo1, fill=gravity)) +
  geom_density(alpha=0.4) + theme_bw()
ggsave("PBMC_unstim_SenMayo.pdf", width = 6, height = 5)
ggplot(pbmc_module@meta.data) + 
  geom_density_ridges(aes(x=SenMayo1, y= predicted.celltype.l2, fill=gravity), alpha=0.5) + theme_bw()
ggsave("PBMC_unstim_SenMayo_celltypes.pdf", width = 6, height = 6)
t.test(filter(pbmc_module@meta.data,gravity=="1G")$SenMayo1,
       filter(pbmc_module@meta.data,gravity=="uG")$SenMayo1)
  # make a dataframe - dat for violin and boxplot
c1<-c(rep("1G",length(filter(pbmc_module@meta.data,gravity=="1G")$SenMayo1)),
      rep("uG",length(filter(pbmc_module@meta.data,gravity=="uG")$SenMayo1)))
c2<-c(filter(pbmc_module@meta.data,gravity=="1G")$SenMayo1,filter(pbmc_module@meta.data,gravity=="uG")$SenMayo1)
dat<-data.frame(Group=c1,Score=c2)
  # violin plot
ggplot(dat,aes(Group,Score)) + 
  geom_violin(aes(fill = Group)) +
  theme_bw() + scale_fill_brewer(palette = 'Set2',direction = -1) +
  labs(title="SenMayo") +
  geom_boxplot(width=0.1,outlier.shape = NA) +
  ggpubr::stat_compare_means(method = "t.test", 
                             comparisons=list(c("1G","uG")), 
                             label = 'p.signif')
ggsave("uGvs1G_unstim_SenMayo_violin.pdf", width = 5, height = 7)
  # boxplot
ggboxplot(dat, x="Group",y="Score",
          color = "Group", palette = "jco",outlier.shape = NA) + coord_cartesian(ylim = c(-0.1,0.15)) +
  stat_compare_means(method = "wilcox.test", label = "p.signif",label.y=0.11,
                     comparisons=list(c("1G","uG")))
ggsave("uGvs1G_unstim_SenMayo_box.pdf", width = 6, height = 5)


## Reclustering for comparing with Twins study with 4 cell types
# check markers CD4, CD8, CD19, and LD
DefaultAssay(Sample.combined)
Idents(Sample.combined)
unique(Sample.combined$predicted.celltype.l2)
VlnPlot(Sample.combined, features = c("CD4","CD8A"))
VlnPlot(Sample.combined, features = c("CD19"))
Sample.combined$CellType4Twins <- case_when(Sample.combined$predicted.celltype.l2 %in%
                                    c("CD4 Naive","CD4 TCM","CD4 TEM","Treg","CD4 CTL","CD4 Proliferating") ~ "CD4",
                                  Sample.combined$predicted.celltype.l2 %in%
                                    c("MAIT","CD8 TEM","CD8 Naive","CD8 TCM") ~ "CD8",
                                  Sample.combined$predicted.celltype.l2 %in%
                                    c("B memory","B intermediate","B naive") ~ "CD19",
                                  Sample.combined$predicted.celltype.l2 %in%
                                    c("CD14 Mono","CD16 Mono","cDC1","cDC2","pDC","ASDC") ~ "LD",
                                  TRUE ~ as.character(Sample.combined$predicted.celltype.l2))

Sample.combined$CellType4Twins_treat <- paste0(Sample.combined$CellType4Twins,"_",Sample.combined$treatment)

Idents(Sample.combined) <- "CellType4Twins_treat"
plan("multisession", workers = 8)

CD4_uGvs1G_unstim <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="CD4_uG_unstimulated", ident.2 ="CD4_1G_unstimulated",
                          min.pct=0.005, logfc.threshold = 0.1,
                          test.use = "MAST")

CD8_uGvs1G_unstim <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="CD8_uG_unstimulated", ident.2 ="CD8_1G_unstimulated",
                                 min.pct=0.005, logfc.threshold = 0.1,
                                 test.use = "MAST")

CD19_uGvs1G_unstim <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="CD19_uG_unstimulated", ident.2 ="CD19_1G_unstimulated",
                                 min.pct=0.005, logfc.threshold = 0.1,
                                 test.use = "MAST")

LD_uGvs1G_unstim <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="LD_uG_unstimulated", ident.2 ="LD_1G_unstimulated",
                                  min.pct=0.005, logfc.threshold = 0.1,
                                  test.use = "MAST")

# save the marker lists for CD4, CD8, CD19, and LD
write.csv(CD4_uGvs1G_unstim, "CD4_uGvs1G_unstim.csv")
write.csv(CD8_uGvs1G_unstim, "CD8_uGvs1G_unstim.csv")
write.csv(CD19_uGvs1G_unstim, "CD19_uGvs1G_unstim.csv")
write.csv(LD_uGvs1G_unstim, "LD_uGvs1G_unstim.csv")


### --- collect normalized counts in SC pseudo-bulk for comparision with real bulk RNAseq --- ##
# change/check analysis slots
Idents(Sample.combined) <- "treatment"
DefaultAssay(Sample.combined)
# uG
cell.uG.unstim <- Sample.combined@meta.data %>% filter(treatment=="uG_unstimulated") %>% row.names()
SCT.uG.unstim <- Sample.combined@assays[["SCT"]]@data[,cell.uG.unstim]
#write.csv(SCT.uG.unstim,"SCT_uG_unstim.csv")

SCT.uG.unstim.mean <- data.frame(gene_name=row.names(SCT.uG.unstim),SC_uG_means=rowSums(SCT.uG.unstim)/2)

bulk.full <- read.csv("~/scRNAseq_analysis/host_counts_DEG.csv")
bulk.test.uG <- bulk.full %>% dplyr::select(gene_name=hybrid_name,uG_0720.norm,uG_0807.norm,uG_0810.norm)
bulk.test.uG$bulk_mean <- rowMeans(bulk.test.uG[2:4])

SC.bulk.test.uG <- inner_join(SCT.uG.unstim.mean,bulk.test.uG, by="gene_name")

# 1G
cell.1G.unstim <- Sample.combined@meta.data %>% filter(treatment=="1G_unstimulated") %>% row.names()
SCT.1G.unstim <- Sample.combined@assays[["SCT"]]@data[,cell.1G.unstim]
#write.csv(SCT.1G.unstim,"SCT_1G_unstim.csv")

SCT.1G.unstim.mean <- data.frame(gene_name=row.names(SCT.1G.unstim),SC_1G_means=rowSums(SCT.1G.unstim)/2)

bulk.full <- read.csv("~/scRNAseq_analysis/host_counts_DEG.csv")
bulk.test.1G <- bulk.full %>% dplyr::select(gene_name=hybrid_name,h1G_0720.norm,h1G_0807.norm,h1G_0810.norm)
bulk.test.1G$bulk_mean <- rowMeans(bulk.test.1G[2:4])

SC.bulk.test.1G <- inner_join(SCT.1G.unstim.mean,bulk.test.1G, by="gene_name")

# # check distributions
# ggplot(SC.bulk.test.1G, aes(x=SC_1G_means)) +
#   geom_density()
# ggplot(SC.bulk.test.1G, aes(x=bulk_mean)) +
#   geom_density()

# log10 the data for correlation
## 1G
SC.bulk.test.1G$log_SC_1G_means <- log10(SC.bulk.test.1G$SC_1G_means)
SC.bulk.test.1G$log_bulk_1G_means <- log10(SC.bulk.test.1G$bulk_mean)
ggscatter(SC.bulk.test.1G, x = "log_SC_1G_means", y = "log_bulk_1G_means", 
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "log10 of SC 1G normalized counts", ylab = "log10 of Bulk 1G normalized counts", title = "SC and Bulk 1G",
          size = 0.5,color = brewer.pal(9, 'Paired')[1]) +
  geom_smooth(method="lm",color = brewer.pal(9, 'Paired')[2])
ggsave("SCXBulk_1G_spearman_log.pdf",width = 5, height = 4)

## uG
SC.bulk.test.uG$log_SC_uG_means <- log10(SC.bulk.test.uG$SC_uG_means)
SC.bulk.test.uG$log_bulk_uG_means <- log10(SC.bulk.test.uG$bulk_mean)
ggscatter(SC.bulk.test.uG, x = "log_SC_uG_means", y = "log_bulk_uG_means", 
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "log10 of SC uG normalized counts", ylab = "log10 of Bulk uG normalized counts", title = "SC and Bulk uG",
          size = 0.5,color = brewer.pal(9, 'Paired')[1]) +
  geom_smooth(method="lm",color = brewer.pal(9, 'Paired')[2])
ggsave("SCXBulk_uG_spearman_log.pdf",width = 5, height = 4)
