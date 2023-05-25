# Include code used in supplementary Figure 12-13
### Setup the Seurat objects ###
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(future)
library(stringr)
# change the current plan to access parallelization
plan("multisession", workers = 16)

setwd("~/scRNAseq_analysis/PBMC_microgravity")

Idents(Sample.combined)<-"treatment"
# stim_vs_unstim_1G overall
stim_vs_unstim_1G <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="1G_stimulated", ident.2 ="1G_unstimulated",
                                min.pct=0.005, logfc.threshold = 0.1,
                                test.use = "MAST")
stim_vs_unstim_1G <- stim_vs_unstim_1G[order(stim_vs_unstim_1G$avg_log2FC, decreasing = T),] # sort by log2FC
write.csv(stim_vs_unstim_1G, "DEGs_overall_stim_vs_unstim_1G.csv")

stim_vs_unstim_1G.all <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="1G_stimulated", ident.2 ="1G_unstimulated",
                                 min.pct=0, logfc.threshold = 0,
                                 test.use = "MAST")
stim_vs_unstim_1G.all <- stim_vs_unstim_1G.all[order(stim_vs_unstim_1G.all$avg_log2FC, decreasing = T),] # sort by log2FC
write.csv(stim_vs_unstim_1G.all, "DEGs_overall_stim_vs_unstim_1G_all.csv")

# stim_vs_unstim_uG overall
stim_vs_unstim_uG <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="uG_stimulated", ident.2 ="uG_unstimulated",
                                 min.pct=0.005, logfc.threshold = 0.1,
                                 test.use = "MAST")
stim_vs_unstim_uG <- stim_vs_unstim_uG[order(stim_vs_unstim_uG$avg_log2FC, decreasing = T),] # sort by log2FC
write.csv(stim_vs_unstim_uG, "DEGs_overall_stim_vs_unstim_uG.csv")

stim_vs_unstim_uG.all <- FindMarkers(Sample.combined, assay = "SCT", ident.1 ="uG_stimulated", ident.2 ="uG_unstimulated",
                                     min.pct=0, logfc.threshold = 0,
                                     test.use = "MAST")
stim_vs_unstim_uG.all <- stim_vs_unstim_uG.all[order(stim_vs_unstim_uG.all$avg_log2FC, decreasing = T),] # sort by log2FC
write.csv(stim_vs_unstim_uG.all, "DEGs_overall_stim_vs_unstim_uG_all.csv")


# subset 1G samples
Idents(Sample.combined)<-"treatment"
Sample.combined.1G <- subset(Sample.combined, idents = c("1G_unstimulated","1G_stimulated"))
# update the levels of treatment
Sample.combined.1G$treatment<-factor(Sample.combined.1G$treatment, 
                                         levels=c("1G_unstimulated","1G_stimulated"))
## view&save total umap in unstimulated subset
Idents(Sample.combined.1G) <- "predicted.celltype.l2"
DimPlot(Sample.combined.1G,split.by = "treatment", label =T, repel = T)
ggsave("umap_predicted.celltype.1G.pdf",width=12,height=8)

# subset uG samples
Idents(Sample.combined)<-"treatment"
Sample.combined.uG <- subset(Sample.combined, idents = c("uG_unstimulated","uG_stimulated"))
# update the levels of treatment
Sample.combined.uG$treatment<-factor(Sample.combined.uG$treatment, 
                                     levels=c("uG_unstimulated","uG_stimulated"))
## view&save total umap in unstimulated subset
Idents(Sample.combined.uG) <- "predicted.celltype.l2"
DimPlot(Sample.combined.uG,split.by = "treatment", label =T, repel = T)
ggsave("umap_predicted.celltype.uG.pdf",width=12,height=8)

# cell proportion changes - 1G
sub.prop1G<-sub.prop.all[sub.prop.all$sample %in% c("1G_unstimulated","1G_stimulated"),]
# cell type ratio change: unstimulated
ggplot(sub.prop1G, aes(x = Var1, y = Freq, fill = sample)) +
  geom_bar(stat = "identity",position = "dodge") +
  ggtitle("Cell type frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(x="Cell types", y="Proportion", title="Cell type frequency") +
  guides(fill=guide_legend(title='Groups'))
ggsave("Cell_type_frequency_1G.pdf",width=10,height=6)

# cell proportion changes - uG
sub.propuG<-sub.prop.all[sub.prop.all$sample %in% c("uG_unstimulated","uG_stimulated"),]
# cell type ratio change: unstimulated
ggplot(sub.propuG, aes(x = Var1, y = Freq, fill = sample)) +
  geom_bar(stat = "identity",position = "dodge") +
  ggtitle("Cell type frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust=1)) +
  labs(x="Cell types", y="Proportion", title="Cell type frequency") +
  guides(fill=guide_legend(title='Groups'))
ggsave("Cell_type_frequency_uG.pdf",width=10,height=6)

# cell type ratio change: log2FC comparision unstim; remove celltype ratio <1%
# Stimulated vs Unstimulated 1G
sub.prop1G.1<-sub.prop1G[sub.prop1G$sample=="1G_unstimulated",]
sub.prop1G.1<-sub.prop1G.1[sub.prop1G.1$Freq>0.01,]
sub.prop1G.2<-sub.prop1G[sub.prop1G$sample=="1G_stimulated",]
sub.prop1G.2<-sub.prop1G.2[sub.prop1G.2$Freq>0.01,]
sub.prop.comm<-intersect(sub.prop1G.1$Var1,sub.prop1G.2$Var1)
sub.prop1G.1<-sub.prop1G.1[sub.prop1G.1$Var1 %in% sub.prop.comm,]
sub.prop1G.2<-sub.prop1G.2[sub.prop1G.2$Var1 %in% sub.prop.comm,]
sub.prop2vs1<-data.frame(Cell_type=sub.prop1G.1$Var1,
                         group="Stimulated vs Unstimulated 1G",
                         Fold_change=sub.prop1G.2$Freq/sub.prop1G.1$Freq)
sub.prop2vs1$log2FC<-log2(sub.prop2vs1$Fold_change)
ggplot(sub.prop2vs1, aes(Cell_type,log2FC))+
  geom_bar(position="dodge", stat="identity",fill="#999999") +
  ggtitle("Stimulated vs Unstimulated in 1G") +
  geom_text(aes(label=round(log2FC,2)),size=2.5,
            position=position_dodge(0.9), vjust=-0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust =0.5),
        axis.text.x = element_text(angle=45, hjust=1))
ggsave("CellType_ratio_Stimulated_vs_Unstimulated_1G.pdf",width=6,height=4)

# Stimulated vs Unstimulated uG
sub.propuG.1<-sub.propuG[sub.propuG$sample=="uG_unstimulated",]
sub.propuG.1<-sub.propuG.1[sub.propuG.1$Freq>0.01,]
sub.propuG.2<-sub.propuG[sub.propuG$sample=="uG_stimulated",]
sub.propuG.2<-sub.propuG.2[sub.propuG.2$Freq>0.01,]
sub.prop.comm<-intersect(sub.propuG.1$Var1,sub.propuG.2$Var1)
sub.propuG.1<-sub.propuG.1[sub.propuG.1$Var1 %in% sub.prop.comm,]
sub.propuG.2<-sub.propuG.2[sub.propuG.2$Var1 %in% sub.prop.comm,]
sub.prop2vs1<-data.frame(Cell_type=sub.propuG.1$Var1,
                         group="Stimulated vs Unstimulated uG",
                         Fold_change=sub.propuG.2$Freq/sub.propuG.1$Freq)
sub.prop2vs1$log2FC<-log2(sub.prop2vs1$Fold_change)
ggplot(sub.prop2vs1, aes(Cell_type,log2FC))+
  geom_bar(position="dodge", stat="identity",fill="#999999") +
  ggtitle("Stimulated vs Unstimulated in uG") +
  geom_text(aes(label=round(log2FC,2)),size=2.5,
            position=position_dodge(0.9), vjust=-0.5) +
  theme_bw() +
  theme(plot.title = element_text(hjust =0.5),
        axis.text.x = element_text(angle=45, hjust=1))
ggsave("CellType_ratio_Stimulated_vs_Unstimulated_uG.pdf",width=6,height=4)

# volcano plot for overall Stimulated vs Unstimulated in 1G
library(EnhancedVolcano)
stim_vs_unstim_1G.all.notaxid<-stim_vs_unstim_1G.all[grep("taxid|MT-",row.names(stim_vs_unstim_1G.all),invert = T),]

EnhancedVolcano(stim_vs_unstim_1G.all.notaxid,
                lab = rownames(stim_vs_unstim_1G.all.notaxid),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                ylab = '-log10(adj.p)',
                legendLabels = c("NS", expression(Log[2] ~ FC), "adj.p-value", expression(adj.p - value ~ and
                                                                                      ~ log[2] ~ FC)),
                xlim =c(-2,3),
                FCcutoff = 0.25,
                labSize = 5, 
                drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                title = "Stimulated vs Unstimulated",
                subtitle = bquote(italic("1G PBMC")))
ggsave("Stimulated_vs_Unstimulated_1G_Volcano.pdf",width=10,height=10)

# volcano plot for overall Stimulated vs Unstimulated in uG
library(EnhancedVolcano)
stim_vs_unstim_uG.all.notaxid<-stim_vs_unstim_uG.all[grep("taxid|MT-",row.names(stim_vs_unstim_uG.all),invert = T),]

EnhancedVolcano(stim_vs_unstim_uG.all.notaxid,
                lab = rownames(stim_vs_unstim_uG.all.notaxid),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                ylab = '-log10(adj.p)',
                legendLabels = c("NS", expression(Log[2] ~ FC), "adj.p-value", expression(adj.p - value ~ and
                                                                                          ~ log[2] ~ FC)),
                xlim =c(-2,3),
                FCcutoff = 0.25,
                labSize = 5, 
                drawConnectors = T, arrowheads=F, min.segment.length=0.3,
                title = "Stimulated vs Unstimulated",
                subtitle = bquote(italic("uG PBMC")))
ggsave("Stimulated_vs_Unstimulated_uG_Volcano.pdf",width=10,height=10)

# save the DEGs for each cell types
# ------ find DEG for all cell types stim vs unstim in 1G ------ #
# with default filter min.pct=0.1 and logfc=0.25
Idents(Sample.combined)<-"celltype.treatment"
plan("multisession", workers = 16)
options(future.globals.maxSize= +Inf) # increase maxSize of future function
cell4DEG<-unique(Sample.combined.unstim$predicted.celltype.l2)
for (d in cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating","gdT","Treg","CD4 CTL")]){
  DEG.1G <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_1G_stimulated"), ident.2 =paste0(d,"_1G_unstimulated"))
  write.csv(DEG.1G, paste0(d,"_stim_vs_unstim_1G.csv"))
}

#### Dot plot:
#### draw dot plot for markers by Log2FC & p-value values
# Stim_vs_Unstim_1G
setwd("Stim_vs_Unstim_1G")
features4plot <- c(head(row.names(stim_vs_unstim_1G),25),tail(row.names(stim_vs_unstim_1G),25))
celltype4dot.1G <- cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating","gdT","Treg","CD4 CTL")]
celltype4dot.1G <- append("Overall",celltype4dot.1G)
lsmarker.dot<-data.frame()
for (d in celltype4dot.1G){
  markers4dot <- read.csv(paste0(d,"_stim_vs_unstim_1G.csv"))
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
  colorspace::scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-3,3), oob = scales::squish, name = 'Log2FC',
                                                na.value="transparent")
ggsave("Dot_Markers_stim_vs_unstim_1G.pdf", width = 11, height = 6,limitsize = FALSE)

# ------ find DEG for all cell types stim vs unstim in uG ------ #
# with default filter min.pct=0.1 and logfc=0.25
Idents(Sample.combined)<-"celltype.treatment"
plan("multisession", workers = 16)
options(future.globals.maxSize= +Inf) # increase maxSize of future function
cell4DEG<-unique(Sample.combined.unstim$predicted.celltype.l2)
for (d in cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating","gdT","Treg","CD4 CTL")]){
  DEG.1G <- FindMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_uG_stimulated"), ident.2 =paste0(d,"_uG_unstimulated"))
  write.csv(DEG.1G, paste0(d,"_stim_vs_unstim_uG.csv"))
}

#### Dot plot:
#### draw dot plot for markers by Log2FC & p-value values
# Stim_vs_Unstim_uG
setwd("Stim_vs_Unstim_uG")
features4plot.uG <- c(head(row.names(stim_vs_unstim_uG.all.notaxid),25),tail(row.names(stim_vs_unstim_uG.all.notaxid),25))
features4plot.uG <- features4plot # select the same gene as 1G
celltype4dot.uG <- cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating","gdT","Treg","CD4 CTL")]
celltype4dot.uG <- append("Overall",celltype4dot.uG)
lsmarker.dot.uG<-data.frame()
for (d in celltype4dot.uG){
  markers4dot <- read.csv(paste0(d,"_stim_vs_unstim_uG.csv"))
  if (nrow(markers4dot) != 0){
    genes<-filter(markers4dot,X %in% features4plot.uG)
    genes<-genes[order(genes$avg_log2FC,decreasing = T),]
    genes$Cell_type<-d
    lsmarker.dot.uG<-rbind(lsmarker.dot.uG,genes)
  }
}
# replace the 0 p-values with lowest number
lsmarker.dot.uG$p_val_adj[lsmarker.dot.uG$p_val_adj==0] <- sort(unique(lsmarker.dot.uG$p_val_adj))[2]

gs<-lsmarker.dot.uG$X %>% unique()
lsmarker.dot.uG$Cell_type<-factor(lsmarker.dot.uG$Cell_type, levels = unique(lsmarker.dot.uG$Cell_type)) # order the x axis

lsmarker.dot.uG %>% filter(X %in% gs) %>% filter(p_val_adj <0.05) %>%
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
  colorspace::scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-3,3), oob = scales::squish, name = 'Log2FC',
                                                na.value="transparent")
ggsave("Dot_Markers_stim_vs_unstim_uG.pdf", width = 11, height = 6,limitsize = FALSE)

# substract: uG -1G
#lsmarker.dot.1G <- lsmarker.dot
#lsmarker.dot.uG_1G <- merge()

stim_vs_unstim_1G.sig <- stim_vs_unstim_1G %>% filter(p_val_adj<0.05)
stim_vs_unstim_uG.sig <- stim_vs_unstim_uG %>% filter(p_val<0.05)

stim_vs_unstim_1GXuG <- merge(stim_vs_unstim_1G.sig,stim_vs_unstim_uG.sig, by="row.names")

stim_vs_unstim_1GXuG <- stim_vs_unstim_1GXuG[order(stim_vs_unstim_1GXuG$avg_log2FC.x, decreasing = T),] # sort by log2FC

# # to determine the top25 based on overall group and the genes of interest for subtract comparison
# gene2check <- c("IFNA2","IFNG","IL1A","IL1B","IL1RA","IL2",
#                 "IL6","IL7","CXCL8","IL17A","IL25","IL17F","CCL2","CCL7","CCL4",
#                 "GBP1","GBP4","GBP5")
# gene2check <- union(head(row.names(stim_vs_unstim_uG.all.notaxid),25),gene2check) # gene in dot plot + genes of interest (ELISA)
# #features4plot.delta <- union(head(row.names(stim_vs_unstim_1GXuG),25),gene2check)
features4plot.delta <- head(row.names(stim_vs_unstim_1G.all.notaxid),50) # unbiased select top50
#features4plot.bot.delta <- tail(row.names(stim_vs_unstim_1GXuG),50)
celltype4dot.delta <- cell4DEG[!cell4DEG %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating","gdT","Treg","CD4 CTL")]
celltype4dot.delta <- append("Overall",celltype4dot.delta)
# 1G DEGs
setwd("../Stim_vs_Unstim_1G")
lsmarker.1G<-list()
for (d in celltype4dot.delta){
  markers4dot <- read.csv(paste0(d,"_stim_vs_unstim_1G.csv"))
  if (nrow(markers4dot) != 0){
    genes<-filter(markers4dot,X %in% features4plot.delta) %>% select(X,avg_log2FC)
    names(genes)[2]<-d
    lsmarker.1G[[d]] <-genes
  }
}
library(purrr)
lsmarker.1G.df <- lsmarker.1G %>% reduce(full_join)

# uG DEGs
setwd("../Stim_vs_Unstim_uG")
lsmarker.uG<-list()
for (d in celltype4dot.delta){
  markers4dot <- read.csv(paste0(d,"_stim_vs_unstim_uG.csv"))
  if (nrow(markers4dot) != 0){
    genes<-filter(markers4dot,X %in% features4plot.delta) %>% select(X,avg_log2FC)
    names(genes)[2]<-d
    lsmarker.uG[[d]] <-genes
  }
}
library(purrr)
lsmarker.uG.df <- lsmarker.uG %>% reduce(full_join)
library(tibble)
lsmarker.1G.df <- column_to_rownames(lsmarker.1G.df,"X")
lsmarker.uG.df <- column_to_rownames(lsmarker.uG.df,"X")

# replace all NA in a dataframe to 0
lsmarker.1G.df[is.na(lsmarker.1G.df)] <- 0
lsmarker.uG.df[is.na(lsmarker.uG.df)] <- 0

lsmarker.1G.df<-lsmarker.1G.df[row.names(lsmarker.uG.df),] # same number of rows
# subtract log2FC values between uG and 1G
lsmarker.delta.uG_1G.df <- lsmarker.uG.df - lsmarker.1G.df

library(ComplexHeatmap)
library(tibble)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(file="Delta_DEGs_stim_vs_unstim_uG_1G_top50_cluster.pdf", width=6, height = 10)
Heatmap(lsmarker.delta.uG_1G.df, col=col_fun,
        cluster_columns = T,
        heatmap_legend_param = list(title = "delta\nuG-1G"))
dev.off()


# IPA pathway delta comparision
library(readxl)
library(tidyr)
library(colorspace)
# read the Excel exported from IPA comparison
IPA_Z_1G <- read_excel("/opt/home/buckcenter.org/fwu/scRNAseq_analysis/PBMC_microgravity/stim_vs_unstim_1G_IPA_Pathways_comparison.xls",skip=1)
IPA_Z_uG <- read_excel("/opt/home/buckcenter.org/fwu/scRNAseq_analysis/PBMC_microgravity/stim_vs_unstim_uG_IPA_Pathways_comparison.xls",skip=1)

# select top50 pathways by z-scores
IPA_Z_1G <- IPA_Z_1G[1:50,]
IPA_Z_uG <- IPA_Z_uG %>% filter(`Canonical Pathways` %in% IPA_Z_1G$`Canonical Pathways`)
# replace Greek symbols
library(stringi)
Greek <- c("α","β","γ","κ","θ")
English <- c("a","b","r","k","th")
IPA_Z_1G$`Canonical Pathways` <- stri_replace_all_regex(IPA_Z_1G$`Canonical Pathways`,
                                                     pattern=Greek,
                                                     replacement = English,
                                                     vectorize=F)
IPA_Z_uG$`Canonical Pathways` <- stri_replace_all_regex(IPA_Z_uG$`Canonical Pathways`,
                                                        pattern=Greek,
                                                        replacement = English,
                                                        vectorize=F)
# Order data frame rows in IPA_Z_uG according to the order in IPA_Z_1G
IPA_Z_uG<-IPA_Z_uG[match(IPA_Z_1G$'Canonical Pathways',IPA_Z_uG$'Canonical Pathways'),]

# trim the name of cell types
colnames(IPA_Z_1G) <- gsub("_stim_vs_unstim.*","",colnames(IPA_Z_1G))
colnames(IPA_Z_uG) <- gsub("_stim_vs_unstim.*","",colnames(IPA_Z_uG))

IPA_Z_1G[IPA_Z_1G=="N/A"] <- NA
IPA_Z_uG[IPA_Z_uG=="N/A"] <- NA

# replace the different pathway name, one of pathway is missing/NA in IPA_Z_uG
IPA_Z_uG$`Canonical Pathways`[is.na(IPA_Z_uG$`Canonical Pathways`)] <- setdiff(IPA_Z_1G$`Canonical Pathways`,IPA_Z_uG$`Canonical Pathways`)

IPA_Z_1G <- column_to_rownames(IPA_Z_1G,'Canonical Pathways')
IPA_Z_uG <- column_to_rownames(IPA_Z_uG,'Canonical Pathways')

# convert the whole dataframe to numeric
IPA_Z_1G <- IPA_Z_1G %>% mutate_all(as.numeric)
IPA_Z_uG <- IPA_Z_uG %>% mutate_all(as.numeric)

# replace all NA in a dataframe to 0
IPA_Z_1G[is.na(IPA_Z_1G)] <- 0
IPA_Z_uG[is.na(IPA_Z_uG)] <- 0

# subtract the IPA z-score values: uG - 1G
IPA_Z_delta <- IPA_Z_uG - IPA_Z_1G

# delta IPA heatmap
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

pdf(file="Delta_IPA_stim_vs_unstim_uG_1G.pdf", width=11, height = 10)
Heatmap(IPA_Z_delta, 
        col=col_fun,
        cluster_columns = F,
        cluster_rows = T,
        row_names_max_width = max_text_width(rownames(IPA_Z_delta)),
        heatmap_legend_param = list(title = "delta\nuG-1G"))
dev.off()

# export to Tab-delimited tabular input format (.txt)
Sample.combined@assays[["RNA"]]@counts
Sample.combined$predicted.celltype.l2

df1 <- as.data.frame(Sample.combined@assays[["RNA"]]@counts)
df2 <- Sample.combined@meta.data %>% select(predicted.celltype.l2)
mapping <- setNames(as.vector(df2[,1]), rownames(df2))
# Change the column names of df1 based on the mapping
colnames(df1) <- sapply(colnames(df1), function(x) ifelse(x %in% names(mapping), mapping[x], x))

#devtools::install_github("RRHO2/RRHO2", build_opts = c("--no-resave-data", "--no-manual"))
library(RRHO2)
# Create "gene" lists:
stim_vs_unstim_1G.all.rank <- stim_vs_unstim_1G.all %>% 
  select(avg_log2FC,p_val_adj) %>%
  mutate(p_val_adj=ifelse(p_val_adj==0,1e-305,p_val_adj)) %>%
  rownames_to_column("gene_symbol")
stim_vs_unstim_1G.all.rank$weight <- stim_vs_unstim_1G.all.rank$avg_log2FC*-log10((stim_vs_unstim_1G.all.rank$p_val_adj))
stim_vs_unstim_1G.all.rank$rank <- order(stim_vs_unstim_1G.all.rank$weight,decreasing=T)
stim_vs_unstim_1G.all.rank <- stim_vs_unstim_1G.all.rank %>% select(gene_symbol,weight)

stim_vs_unstim_uG.all.rank <- stim_vs_unstim_uG.all %>% 
  select(avg_log2FC,p_val_adj) %>%
  mutate(p_val_adj=ifelse(p_val_adj==0,1e-305,p_val_adj)) %>%
  rownames_to_column("gene_symbol")
stim_vs_unstim_uG.all.rank$weight <- stim_vs_unstim_uG.all.rank$avg_log2FC*-log10((stim_vs_unstim_uG.all.rank$p_val_adj))
stim_vs_unstim_uG.all.rank$rank <- order(stim_vs_unstim_uG.all.rank$weight,decreasing=T)
stim_vs_unstim_uG.all.rank <- stim_vs_unstim_uG.all.rank %>% select(gene_symbol,weight)

# Compute overlap and significance
RRHO_obj <-  RRHO2_initialize(stim_vs_unstim_1G.all.rank, stim_vs_unstim_uG.all.rank, labels = c("Stim_vs_unstim in 1G", "Stim_vs_unstim in uG"), log10.ind=TRUE)

pdf(file="RRHO2_stim_vs_unstim_1GuG.pdf", width=7, height = 6)
RRHO2_heatmap(RRHO_obj)
dev.off()

pdf(file="RRHO2_vennDiagra_DD_stim_vs_unstim_1GuG.pdf", width=6, height = 5)
RRHO2_vennDiagram(RRHO_obj,type="dd")
dev.off()

pdf(file="RRHO2_vennDiagra_UU_stim_vs_unstim_1GuG.pdf", width=6, height = 5)
RRHO2_vennDiagram(RRHO_obj,type="uu")
dev.off()

# extract the uu gene expressions
uu<-RRHO_obj[["genelist_uu"]][["gene_list_overlap_uu"]]
# extract the dd gene expressions
dd<-RRHO_obj[["genelist_dd"]][["gene_list_overlap_dd"]]

# uu in sn
uu.1G<-select(stim_vs_unstim_1G.all,c(avg_log2FC,p_val_adj))[uu,]
colnames(uu.1G)<-paste0(colnames(uu.1G),"_1G")
# uu in GeoMx
uu.uG<-select(stim_vs_unstim_uG.all,c(avg_log2FC,p_val_adj))[uu,]
colnames(uu.uG)<-paste0(colnames(uu.uG),"_uG")

uu.1GuG <- cbind(uu.1G,uu.uG)
write.csv(uu.1GuG,"RRHO_stim_vs_unstim_uu_1GuG.csv")

# dd in sn
dd.1G<-select(stim_vs_unstim_1G.all,c(avg_log2FC,p_val_adj))[dd,]
colnames(dd.1G)<-paste0(colnames(dd.1G),"_1G")
# dd in GeoMx
dd.uG<-select(stim_vs_unstim_uG.all,c(avg_log2FC,p_val_adj))[dd,]
colnames(dd.uG)<-paste0(colnames(dd.uG),"_uG")

dd.1GuG <- cbind(dd.1G,dd.uG)
write.csv(dd.1GuG,"RRHO_stim_vs_unstim_dd_1GuG.csv")


### conserved DEGs
Sample.combined@meta.data <- tidyr::separate(Sample.combined@meta.data,treatment,sep="_",remove=F,
         into=c("gravity","stimulus"))
Idents(Sample.combined) <- "gravity"
# find the genes that conserved in between uG and 1G conditions, irrespective of stimulus
conserved_uG <- FindConservedMarkers(Sample.combined, assay = "SCT", ident.1 ="uG", ident.2 ="1G",
                                       grouping.var = "stimulus",
                                       min.pct=0.005, logfc.threshold = 0.1,
                                       test.use = "MAST")
write.csv(conserved_uG,"Overall_conserved_uGvs1G.csv.csv")


# ------ find conserved DEG for all cell types ------ #
Sample.combined$cell_gravity <- paste0(Sample.combined$predicted.celltype.l2,"_",Sample.combined$gravity)
# with default filter min.pct=0.1 and logfc=0.25
Idents(Sample.combined)<-"cell_gravity"
plan("multisession", workers = 4)
options(future.globals.maxSize= +Inf) # increase maxSize of future function
cell4DEG.conserved<-unique(Sample.combined$predicted.celltype.l2)
sort(table(Sample.combined$cell_gravity)) # check the cell number <3
cell4DEG.conserved<-cell4DEG.conserved[!cell4DEG.conserved %in% c("Unknow","ASDC","cDC1","CD4 Proliferating","ILC","NK Proliferating")]
for (d in cell4DEG.conserved){
  conserved_DEG_cell <- FindConservedMarkers(Sample.combined, assay = "SCT", ident.1 =paste0(d,"_uG"), ident.2 =paste0(d,"_1G"),
                        grouping.var = "stimulus",
                        min.pct=0.005, logfc.threshold = 0.1,
                        test.use = "MAST")
  write.csv(conserved_DEG_cell, paste0("Conserved_DEGs/",d,"_conserved_uGvs1G.csv"))
}

# filter the conserved gene max_pval <0.05
conserved_uG.filter <- conserved_uG %>% 
  filter(max_pval<0.05) %>%
  filter((stimulated_avg_log2FC>0 & unstimulated_avg_log2FC>0)|(stimulated_avg_log2FC<0 & unstimulated_avg_log2FC<0))
# calculate the top genes by sum of log2FC
conserved_uG.filter$sum_log2FC <- abs(conserved_uG.filter$stimulated_avg_log2FC+conserved_uG.filter$unstimulated_avg_log2FC)
# sort by sum_log2FC
conserved_uG.filter<-conserved_uG.filter[order(conserved_uG.filter$sum_log2FC,decreasing =T),]
# select top50 genes as targets for plot
features4plot.conserved <- head(row.names(conserved_uG.filter),50)
celltype4dot.conserved <- append("Overall",cell4DEG.conserved)
lsmarker.dot.c<-data.frame()
for (d in celltype4dot.conserved){
  markers4dot <- read.csv(paste0("Conserved_DEGs/",d,"_conserved_uGvs1G.csv"),row.names=1)
  if (nrow(markers4dot) != 0 & ncol(markers4dot)>5){
    # drop the row contains taxid or MT genes
    markers4dot <- subset(markers4dot, !(grepl("taxid",rownames(markers4dot)) | grepl("MT-",rownames(markers4dot))))
    # select the target conserved genes
    genes<-filter(markers4dot,rownames(markers4dot) %in% features4plot.conserved)
    # select "stimulated_" group
    genes.stim <- genes[1:5]
    colnames(genes.stim) <- gsub("^stimulated_","",colnames(genes.stim))
    genes.stim$stimulus <- "stimulated"
    # select "unstimulated_" group
    genes.unstim <- genes[6:10]
    colnames(genes.unstim) <- gsub("^unstimulated_","",colnames(genes.unstim))
    genes.unstim$stimulus <- "unstimulated"
    # move the row name to a column to avoid the wrong gene names
    genes.stim <- rownames_to_column(genes.stim,"Genes")
    genes.unstim <- rownames_to_column(genes.unstim,"Genes")
    # combine into one dataframe in a longer shape like "pivot_longer"
    genes.rbind <- rbind(genes.unstim,genes.stim)
    # add cell_type column
    genes.rbind$cell_type <- d
    # iterate add for each cell type
    lsmarker.dot.c<-rbind(lsmarker.dot.c,genes.rbind)
  }
}

# replace the 0 p-values with lowest number
lsmarker.dot.c$p_val_adj[lsmarker.dot.c$p_val_adj==0] <- sort(unique(lsmarker.dot.c$p_val_adj))[2]

# add Stimulus_Cell column for plot
lsmarker.dot.c$Stimulus_Cell <- paste0(lsmarker.dot.c$cell_type,"_",lsmarker.dot.c$stimulus)

# re-order the Stimulus_Cell column for the plot
lsmarker.dot.c$Stimulus_Cell <- factor(lsmarker.dot.c$Stimulus_Cell, levels = names(sort(table(lsmarker.dot.c$Stimulus_Cell),decreasing=T)))

# # replace the 0 p-values with .Machine$double.xmin = 2.225074e-308, which is the smallest number that R can represent
# lsmarker.dot.c$p_val_adj[lsmarker.dot.c$p_val_adj==0] <- .Machine$double.xmin

### ..ooOOoo.. dot/bubble plot ..ooOOoo.. ###
# Calculate the y-axis positions for horizontal lines
y_positions <- seq(-1.55, length(unique(lsmarker.dot.c$Stimulus_Cell))+1, by = 2)

lsmarker.dot.c %>%
ggplot(aes(x=Genes, y = Stimulus_Cell, color = avg_log2FC, size = -log10(p_val_adj))) + 
  geom_point() +
  scale_y_discrete(labels=function(y) str_wrap(y, width=70)) +
  ylab('Cell Types') +
  xlab('Genes') +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size=9, angle=45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=10))+
  theme(axis.line  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  #scale_size(breaks=c(0,25,50,100,200)) +
  scale_size_area(
    max_size = 5,
    breaks = c(0,25,50,100),
    labels = c("0","25","50","100+"),
    guide = "legend",
    limits = c(0, 100),
    oob = scales::squish
  )+
  colorspace::scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-1,1), oob = scales::squish, name = 'Log2FC',
                                                na.value="transparent") +
  geom_hline(yintercept = y_positions, linetype = "dashed", color = "gray") +
  ggtitle("Conserved Genes in uG vs 1G") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Dot_Conserved_Markers_uGvs1G_top50.pdf", width = 13, height = 10)


# Save the conserved genes to a single .csv file
lsmarker.dot.c.all<-data.frame()
for (d in celltype4dot.conserved){
  markers4dot <- read.csv(paste0("Conserved_DEGs/",d,"_conserved_uGvs1G.csv"),row.names=1)
  if (nrow(markers4dot) != 0 & ncol(markers4dot)>5){
    # drop the row contains taxid or MT genes
    markers4dot <- subset(markers4dot, !(grepl("taxid",rownames(markers4dot)) | grepl("MT-",rownames(markers4dot))))
    # move the row name to a column to avoid the wrong gene names
    markers4dot <- rownames_to_column(markers4dot,"Genes")
    # add cell_type column
    markers4dot$cell_type <- d
    # iterate add for each cell type
    lsmarker.dot.c.all<-rbind(lsmarker.dot.c.all,markers4dot)
  }
}
write.csv(lsmarker.dot.c.all,"Conserved_DEGs_uGvs1G_all.csv")

# select all the stimulation induced genes from stim vs unstim 1G
stim_vs_unstim_1G.all.notaxid.induced <- stim_vs_unstim_1G.all.notaxid %>% filter(avg_log2FC>0)
write.csv(stim_vs_unstim_1G.all.notaxid.induced,"stim_vs_unstim_1G.all.notaxid.induced.csv")

