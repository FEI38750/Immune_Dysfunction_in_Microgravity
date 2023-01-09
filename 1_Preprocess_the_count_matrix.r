# ------ 1: preprocess the count matrix from 10x Cellranger and MTD pipelines ------ #

### Setup the Seurat objects ###
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(glmGamPoi)
library(future)
library(tibble)
library(dplyr)
library(parallel)
# change the current plan to access parallelization
plan("multisession", workers = 16)

setwd("~/scRNAseq_analysis/PBMC_microgravity")
# import 10X results
  # Female PBMC
F.uG_unstimulated <- Read10X(data.dir = "~/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/Microgravity_uG_un-simulated/outs/filtered_feature_bc_matrix/")
F.uG_stimulated <- Read10X(data.dir = "~/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/Microgravity_uG_simulated/outs/filtered_feature_bc_matrix/")
F.1G_unstimulated <- Read10X(data.dir = "~/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/Control_1G_un-stimulated/outs/filtered_feature_bc_matrix/")
F.1G_stimulated <- Read10X(data.dir = "~/scRNAseq_runs/scRNA_hPBMCmicrogravity_050622/Control_1G_stimulated/outs/filtered_feature_bc_matrix/")
  # Male PBMC
M.uG_unstimulated <- Read10X(data.dir = "~/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/Microgravity_uG_unstimulated/outs/filtered_feature_bc_matrix/")
M.uG_stimulated <- Read10X(data.dir = "~/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/Microgravity_uG_stimulated/outs/filtered_feature_bc_matrix/")
M.1G_unstimulated <- Read10X(data.dir = "~/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/Control_1G_unstimulated/outs/filtered_feature_bc_matrix/")
M.1G_stimulated <- Read10X(data.dir = "~/scRNAseq_runs/scRNA_hPBMCmicrogravity_052522/Control_1G_stimulated/outs/filtered_feature_bc_matrix/")
# make a sample list
Sample.list <- list("F.uG_unstimulated"=F.uG_unstimulated, 
                    "F.uG_stimulated"=F.uG_stimulated, 
                    "F.1G_unstimulated"=F.1G_unstimulated,
                    "F.1G_stimulated"=F.1G_stimulated,
                    "M.uG_unstimulated"=M.uG_unstimulated,
                    "M.uG_stimulated"=M.uG_stimulated,
                    "M.1G_unstimulated"=M.1G_unstimulated,
                    "M.1G_stimulated"=M.1G_stimulated)
# transfer to dataframe for reprocessing
Sample.list <- lapply(Sample.list, as.data.frame)
# clean the names of cell barcodes
for (i in 1:length(Sample.list)){
  colnames(Sample.list[[i]])<-gsub("-1","",colnames(Sample.list[[i]]))
}

# Import MTD result
  # Female PBMC
F.uG_unstimulated.MTD<-read.delim("~/scRNAseq_analysis/MTD_F/output/Microgravity_uG_unstimulated_count_matrix.txt", row.names = 1)
F.uG_stimulated.MTD<-read.delim("~/scRNAseq_analysis/MTD_F/output/Microgravity_uG_stimulated_count_matrix.txt", row.names = 1)
F.1G_unstimulated.MTD<-read.delim("~/scRNAseq_analysis/MTD_F/output/Control_1G_unstimulated_count_matrix.txt", row.names = 1)
F.1G_stimulated.MTD<-read.delim("~/scRNAseq_analysis/MTD_F/output/Control_1G_stimulated_count_matrix.txt", row.names = 1)
  # Male PBMC
M.uG_unstimulated.MTD<-read.delim("~/scRNAseq_analysis/MTD_M/output/Microgravity_uG_unstimulated_count_matrix.txt", row.names = 1)
M.uG_stimulated.MTD<-read.delim("~/scRNAseq_analysis/MTD_M/output/Microgravity_uG_stimulated_count_matrix.txt", row.names = 1)
M.1G_unstimulated.MTD<-read.delim("~/scRNAseq_analysis/MTD_M/output/Control_1G_unstimulated_count_matrix.txt", row.names = 1)
M.1G_stimulated.MTD<-read.delim("~/scRNAseq_analysis/MTD_M/output/Control_1G_stimulated_count_matrix.txt", row.names = 1)
# combine host and microbiome Sc data in dataframe
for (i in names(Sample.list)){
  Sample.list[[i]]<-rbind(Sample.list[[i]],get(paste0(i,".MTD")))
}

# create seuratobject for QC
Sample.list.qc<-list()
for (i in names(Sample.list)){
  Sample.list.qc[[i]] <- CreateSeuratObject(Sample.list[[i]], project=i,min.cells = 3, min.feature = 250)
}
# function for Pre-process Seurat object: QC, PCA and UMAP
SCT_m<-function(l){
  l <- PercentageFeatureSet(l,pattern = "^MT-", col.name = "percent.mt")
  l <- subset(l, subset= nFeature_RNA>250 & percent.mt < 10)
  l <- SCTransform(l, vst.flavor = "v2", method = "glmGamPoi",verbose = FALSE)
  l <- RunPCA(l, verbose = FALSE)
  l <- RunUMAP(l,dims = 1:30, verbose = FALSE)
  l <- FindNeighbors(l, dims = 1:30, verbose = FALSE)
  l <- FindClusters(l, verbose = FALSE)
}

# SCTransform for each dataset independently
Sample.list.qc<-mclapply(Sample.list.qc,SCT_m, mc.cores=8)
# Remove doublet
library(DoubletFinder)
for (i in names(Sample.list.qc)){
  # pK Identification (no ground-truth)
  sweep.res.list_sample <- paramSweep_v3(Sample.list.qc[[i]], PCs = 1:30, sct = T)
  sweep.stats_sample <- summarizeSweep(sweep.res.list_sample, GT = FALSE)
  bcmvn_sample <- find.pK(sweep.stats_sample)
  pK<-as.numeric(as.character(bcmvn_sample$pK))[bcmvn_sample$BCmetric==max(bcmvn_sample$BCmetric)]
  ## Homotypic Doublet Proportion Estimate
  annotations <- Sample.list.qc[[i]]@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(Sample.list.qc[[i]]@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  ## Run DoubletFinder with varying classification stringencies
  Sample.list.qc[[i]] <- doubletFinder_v3(Sample.list.qc[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  re.pAnn<-names(Sample.list.qc[[i]]@meta.data)[(length(Sample.list.qc[[i]]@meta.data)-1)]
  Sample.list.qc[[i]] <- doubletFinder_v3(Sample.list.qc[[i]], PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = re.pAnn, sct = T)
}

# filter samples based on results of QC & Doublet
for (i in names(Sample.list.qc)){
  sample.meta.data<-Sample.list.qc[[i]]@meta.data
  singlet<-row.names(sample.meta.data)[sample.meta.data[length(sample.meta.data)]=="Singlet"]
  Sample.list[[i]]<-Sample.list[[i]][singlet]
}
# create seuratobject for integration
for (i in names(Sample.list)){
  Sample.list[[i]] <- CreateSeuratObject(Sample.list[[i]], project=i)
  Sample.list[[i]] <- SCTransform(Sample.list[[i]], vst.flavor = "v2", method = "glmGamPoi",verbose = F)
}
# Perform integration
features <- SelectIntegrationFeatures(object.list = Sample.list, nfeatures = 3000)
Sample.list <- PrepSCTIntegration(object.list = Sample.list, anchor.features = features)
plan("multisession", workers = 16)
  #options(future.globals.maxSize= +Inf) # increase maxSize of future function
Anchors <- FindIntegrationAnchors(object.list = Sample.list, normalization.method = "SCT",
                                  anchor.features = features)
Sample.combined <- IntegrateData(anchorset = Anchors, normalization.method = "SCT")

# Dimensional reduction
Sample.combined <- RunPCA(Sample.combined, verbose = FALSE)
  # ElbowPlot(Sample.combined,ndims = 50)
Sample.combined <- RunUMAP(Sample.combined, reduction = "pca", dims = 1:30)
# Cluster the cells
Sample.combined <- FindNeighbors(Sample.combined, dims = 1:30)
plan("multisession", workers = 32)
Sample.combined <- FindClusters(Sample.combined, resolution = 0.8)

# Prepare object to run differential expression on SCT assay with multiple models
Sample.combined <- PrepSCTFindMarkers(Sample.combined)
# find markers for every cluster compared to all remaining cells
plan("multisession", workers = 32)
markers.all <- FindAllMarkers(Sample.combined, assay = "SCT", min.pct = 0.1, logfc.threshold = 0.25)

# Cell type annotation by SingleR
library(SingleR)
monaco.ref <- celldex::MonacoImmuneData()
ImmGen.ref <- celldex::ImmGenData()
# set the active assay back to “RNA,” and re-do the normalization
DefaultAssay(Sample.combined) <- "RNA"
Sample.combined <- NormalizeData(Sample.combined)
# convert Seurat object to single cell experiment (SCE) for convenience
sce <- as.SingleCellExperiment(DietSeurat(Sample.combined))
monaco.main <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
monaco.fine <- SingleR(test = sce,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
ImmGen.main <- SingleR(test = sce,assay.type.test = 1,ref = ImmGen.ref,labels = ImmGen.ref$label.main)
ImmGen.fine <- SingleR(test = sce,assay.type.test = 1,ref = ImmGen.ref,labels = ImmGen.ref$label.fine)
# see the summary of general cell type annotations
table(monaco.main$pruned.labels)
table(monaco.fine$pruned.labels)
table(ImmGen.main$pruned.labels)
table(ImmGen.fine$pruned.labels)
# add the annotations to the Seurat object metadata
Sample.combined@meta.data$monaco.main <- monaco.main$pruned.labels
Sample.combined@meta.data$monaco.fine <- monaco.fine$pruned.labels
Sample.combined@meta.data$ImmGen.main <- ImmGen.main$pruned.labels
Sample.combined@meta.data$ImmGen.fine <- ImmGen.fine$pruned.labels

# Import Azimuth annotation results
predictions.azimuth<-data.frame()
predictions.files<-list.files(path = '~/scRNAseq_analysis/PBMC_microgravity',pattern='*_azimuth_pred.tsv',full.names =T)
for (p in predictions.files){
  predictions<-read.delim(p, row.names = 1)
  pred.basename<-basename(p)
  if (pred.basename=="F_Microgravity_uG_un-simulated_azimuth_pred.tsv"){
    rownames(predictions)<-gsub("-[0-9]","_1",rownames(predictions))
  } else if (pred.basename== "F_Microgravity_uG_simulated_azimuth_pred.tsv"){
    rownames(predictions)<-gsub("-[0-9]","_2",rownames(predictions))
  } else if (pred.basename== "F_Control_1G_un-stimulated_azimuth_pred.tsv"){
    rownames(predictions)<-gsub("-[0-9]","_3",rownames(predictions))
  } else if (pred.basename== "F_Control_1G_stimulated_azimuth_pred.tsv"){
    rownames(predictions)<-gsub("-[0-9]","_4",rownames(predictions))
  } else if (pred.basename== "M_Microgravity_uG_unstimulated_azimuth_pred.tsv"){
    rownames(predictions)<-gsub("-[0-9]","_5",rownames(predictions))
  } else if (pred.basename== "M_Microgravity_uG_stimulated_azimuth_pred.tsv"){
    rownames(predictions)<-gsub("-[0-9]","_6",rownames(predictions))
  } else if (pred.basename== "M_Control_1G_unstimulated_azimuth_pred.tsv"){
    rownames(predictions)<-gsub("-[0-9]","_7",rownames(predictions))
  } else if (pred.basename== "M_Control_1G_stimulated_azimuth_pred.tsv"){
    rownames(predictions)<-gsub("-[0-9]","_8",rownames(predictions))
  }
  predictions.azimuth<-rbind(predictions.azimuth,predictions)
}
predictions.azimuth<-predictions.azimuth[row.names(Sample.combined@meta.data),]
Sample.combined <- AddMetaData(Sample.combined, metadata = predictions.azimuth)

# Remove RBC
library(dplyr)
Sample.combined@meta.data <- Sample.combined@meta.data %>% replace(is.na(.), "Unknow") # replace all NA
Idents(Sample.combined)<-"predicted.celltype.l2"
DefaultAssay(Sample.combined) <- "integrated"
Sample.combined <- subset(Sample.combined, idents = c("Eryth","Platelet"), invert = T)
# check new cell numbers
for (i in unique(Sample.combined$orig.ident)){
  print(paste0(i,": ",length(Sample.combined$orig.ident[Sample.combined$orig.ident==i])))
}
DimPlot(Sample.combined, reduction = "umap")

# add more metadata
Sample.combined@meta.data$sex <- case_when(Sample.combined@meta.data$orig.ident ==
                                             "F.uG_unstimulated" ~ "Female",
                                           Sample.combined@meta.data$orig.ident ==
                                             "F.uG_stimulated" ~ "Female",
                                           Sample.combined@meta.data$orig.ident ==
                                             "F.1G_unstimulated" ~ "Female",
                                           Sample.combined@meta.data$orig.ident ==
                                             "F.1G_stimulated" ~ "Female",
                                           Sample.combined@meta.data$orig.ident ==
                                             "M.uG_unstimulated" ~ "Male",
                                           Sample.combined@meta.data$orig.ident ==
                                             "M.uG_stimulated" ~ "Male",
                                           Sample.combined@meta.data$orig.ident ==
                                             "M.1G_unstimulated" ~ "Male",
                                           Sample.combined@meta.data$orig.ident ==
                                             "M.1G_stimulated" ~ "Male")
Sample.combined$celltype.sex <- paste(Sample.combined$predicted.celltype.l2, Sample.combined$sex,
                                           sep = "_")
Sample.combined$treatment <- case_when(Sample.combined@meta.data$orig.ident ==
                                         "F.uG_unstimulated" ~ "uG_unstimulated",
                                       Sample.combined@meta.data$orig.ident ==
                                         "F.uG_stimulated" ~ "uG_stimulated",
                                       Sample.combined@meta.data$orig.ident ==
                                         "F.1G_unstimulated" ~ "1G_unstimulated",
                                       Sample.combined@meta.data$orig.ident ==
                                         "F.1G_stimulated" ~ "1G_stimulated",
                                       Sample.combined@meta.data$orig.ident ==
                                         "M.uG_unstimulated" ~ "uG_unstimulated",
                                       Sample.combined@meta.data$orig.ident ==
                                         "M.uG_stimulated" ~ "uG_stimulated",
                                       Sample.combined@meta.data$orig.ident ==
                                         "M.1G_unstimulated" ~ "1G_unstimulated",
                                       Sample.combined@meta.data$orig.ident ==
                                         "M.1G_stimulated" ~ "1G_stimulated")
Sample.combined$celltype.treatment <- paste(Sample.combined$predicted.celltype.l2, Sample.combined$treatment,
                                      sep = "_")
Sample.combined$celltype.sex.treatment <- paste(Sample.combined$predicted.celltype.l2, Sample.combined$orig.ident,
                                            sep = "_")
