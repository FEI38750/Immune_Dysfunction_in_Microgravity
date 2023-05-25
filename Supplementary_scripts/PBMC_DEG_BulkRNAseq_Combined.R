#args = commandArgs(trailingOnly=TRUE) # Passing arguments to an R script from bash/shell command lines
args2<- "/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Bulk_RNAseq/Combined/samplesheet.csv"
args3<- "9606"
args4<-"/Users/feiwu/Library/CloudStorage/GoogleDrive-feiwu38750@gmail.com/My Drive/RNA-Seq/MTD/HostSpecies.csv"

# make folder DEG to store outputs
setwd("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Bulk_RNAseq/Combined")
#system("mkdir -p DEG") # make new a folder DEG in current working directory

library(DESeq2)
library(tibble)
# Read & preprocess the input file cts before run deseq2 analysis
filename<-"host_counts"


# read the data file
#featureCounts file (eg. host_counts.txt) without quote symbol "\""; mark empty quote as quote=""
cts.F<-read.table("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Bulk_RNAseq/Combined/host_counts_F.txt",
                row.names=1, sep="\t",header=T, quote="")
# drop first 5 columns with information other than counts
cts.F<-cts.F[,-c(1:5)]
cts.M<-read.table("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Bulk_RNAseq/Combined/host_counts_M.txt",
                  row.names=1, sep="\t",header=T, quote="")
# drop first 5 columns with information other than counts
cts.M<-cts.M[,-c(1:5)]

# combine F and M
cts <- merge(cts.F,cts.M, by="row.names")
cts <- column_to_rownames(cts,"Row.names")

# drop rows of NA
cts <- na.omit(cts)

# # drop rows of zero count
# cts<-cts[rowSums(cts[-1])>0,]

# drop rows of sum<5 count
cts<-cts[rowSums(cts[-1])>=5,]

# keep the rows that more than 6 (in total 24) columns >0.
cts <- cts[rowSums(cts > 0) >= 6, ]

# make a folder for outputs
dir.create("Host_DEG", recursive = T)
setwd("Host_DEG")

# read the original samplesheet
coldata0 <- read.csv(args2, header = T, na.strings=c("","NA"))
# extract the contrast information for reference
coldata_vs<- coldata0[c("group1","group2")]
coldata_vs<-coldata_vs[rowSums(is.na(coldata_vs)) == 0,] #remove the NA rows
# read samplesheet as factors (as.is = F) for Deseq2 statistical analysis
coldata_factor <- read.csv(args2, header = T, as.is = F)
coldata_factor[]<-lapply(coldata_factor, factor)
coldata<-coldata_factor[,1:2]
# add paired sample information for paired test
coldata$subject <- factor(c(1,2,3,1,2,3,1,2,3,1,2,3,4,5,6,4,5,6,4,5,6,4,5,6))

# make cts(count matrix) has consistent order with samplesheet/metadata
cts<-cts[coldata0$sample_name]
# load the datastructure to DESeq
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design= ~ subject + group)

# perform the DESeq analysis
dds <- DESeq(dds)


# Data transformation for visualization (normalization included)
#rld<-rlog(dds,blind=F) # regularized log transformation (log2 based)
if (dim(results(dds))[1]  < 1000 || min(colSums(cts !=0)) < 1000){
  vsd<-varianceStabilizingTransformation(dds,blind=F) # vatiance stabilizing transformation
} else {
  vsd<-vst(dds,blind=F)
}
normtrans<-assay(vsd)

# save normalized & transformed data for visualization
write.csv(normtrans,file=paste0(sub(".tsv$|.txt$","",filename),"_normalized_transformed.csv"))

# save normalized (untransformed) data for reference
norm<-counts(dds,normalized=T)
write.csv(norm,file=paste0(sub(".tsv$|.txt$","",filename),"_normalized.csv"))

# merge and add suffixes; normalized and normalized&transformed
merge.nt<-merge(norm,normtrans,by="row.names", suffixes=c(".norm",".normtrans"))

library("biomaRt")
host_sp<-read.csv(args4) # read a list of supported host species
# add annotations; try up to 120 times/20 mins if biomaRt server not response
ensembl <- NULL
attempt <- 0
while ( is.null(ensembl) && attempt <= 120){
  try(
    ensembl <- useMart("ensembl",dataset=host_sp[host_sp$Taxon_ID==args3,2])
  )
  attempt<-attempt+1
  if (attempt > 1){
    print(paste0("Retry to get gene annotations from ensembl: ",attempt," times")) 
  }
  Sys.sleep(10)
}
if (is.null(ensembl)){
  stop("Failed to get gene annotations from ensembl server, please try again later")
}
names(merge.nt)[1]<-"GeneID"
genes <- merge.nt$GeneID
gene_ID <- getBM(filters="ensembl_gene_id",
                 attributes=c("external_gene_name","ensembl_gene_id",
                              "chromosome_name","start_position","end_position",
                              "strand","gene_biotype","description"),
                 values=genes,mart=ensembl)
names(gene_ID)[names(gene_ID)=="external_gene_name"]<-"gene_name" #rename the first column

gene_len<-read.table("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Bulk_RNAseq/Combined/host_counts_M.txt",
                       row.names=1, sep="\t",header=T, quote="")
gene_len<-gene_len.M["Length"]

colnames(gene_len)<-"gene_length"
gene_ID<-merge(gene_ID,gene_len,by.x="ensembl_gene_id", by.y="row.names")

# function of contrast between groups, merge to annotation and count tables, and save to .csv files
comparison<-function(dds,coldata_vs,filename,gene_ID,cts,merge.nt){
  for (i in 1:nrow(coldata_vs)){
  group1<-coldata_vs$group1[i]
  group2<-coldata_vs$group2[i]
  comparison <- results(dds, contrast=c("group", group1, group2))
  comparison_f<-as.data.frame(comparison)
  comparison_f<-comparison_f[order(comparison_f$pvalue),]
  system(paste0("mkdir -p"," ",group1,"_vs_",group2)) #for output file structure
  setwd(paste0(group1,"_vs_",group2))
  write.csv(comparison_f,file=paste0(sub(".tsv$|.txt$","",filename),"_",group1,"_vs_",group2,".csv"))
  setwd("../")
  colnames(comparison_f)<-paste0(colnames(comparison_f),".",group1,"_vs_",group2)
  gene_ID<-merge(gene_ID,comparison_f,by.x="ensembl_gene_id", by.y="row.names") #merge each comparison to annotation
  }
  Anno.merge.all<-merge(cts,gene_ID,by.x="row.names",by.y="ensembl_gene_id")
  Anno.merge.all<-merge(Anno.merge.all,merge.nt,by.x="Row.names",by.y="GeneID")
  names(Anno.merge.all)[1]<-"gene_id" #rename the first column
  Anno.merge.all[Anno.merge.all==""]<-"-" #replace the empty cells to "-"
  Anno.merge.all<-Anno.merge.all[complete.cases(Anno.merge.all[,grep("pvalue|padj",names(Anno.merge.all))]),] # drop pvalue == NA rows; optional
  hybrid_name<-Anno.merge.all$gene_name # add a column of gene_id/name hybrid
  while ("-" %in% hybrid_name){
    hybrid_name[match("-",hybrid_name)] <-
      Anno.merge.all$gene_id[match("-",hybrid_name)]}
  Anno.merge.all<-add_column(Anno.merge.all, hybrid_name, .after = "gene_name")
  write.csv(Anno.merge.all,file=paste0(sub(".tsv$|.txt$","",filename),"_DEG.csv"),row.names=F)
}


# apply the comparison function
comparison(dds, coldata_vs, filename,gene_ID,cts,merge.nt)



