### --- Heatmap of IPA comparison results --- ###
## SC, bulk, Q_treatment, GLDS-420 comparison
library(ComplexHeatmap)
library(tibble)
# load IPA pathway comparison z-score
PW_Z <- read.csv("~/scRNAseq_analysis/PBMC/Comparison_SCXBulk_Pathways_zscore.csv")
PW_Z <- column_to_rownames(PW_Z,"Canonical.Pathways")
PW_Z <- as.matrix(PW_Z)
# load IPA pathway comparison adj.p
PW_P <- read.csv("~/scRNAseq_analysis/PBMC/Comparison_SCXBulk_Pathways_adjP.csv")
PW_P <- column_to_rownames(PW_P,"Canonical.Pathways")
PW_P <- as.matrix(PW_P)

# replace Greek symbols
library(stringi)
Greek <- c("α","β","γ","κ","θ")
English <- c("a","b","r","k","th")
rownames(PW_Z) <- stri_replace_all_regex(rownames(PW_Z),
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)
rownames(PW_P) <- stri_replace_all_regex(rownames(PW_P),
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)

# reorder the p value matrix by the z-score matrix
PW_P <- PW_P[rownames(PW_Z),colnames(PW_Z)]

# draw heatmap
library(circlize)
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))

pdf(file="IPA_heatmap_Q.pdf", width=10, height = 40)
ht <- Heatmap(PW_Z,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "top",
              column_dend_height = unit(18, "mm"),
              heatmap_legend_param = list(title="Activation\nz-score",
                                          at=c(-4,-2,0,2,4),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(PW_Z)),
              row_names_gp = gpar(fontsize = 9.5))
ht<-draw(ht)
dev.off()



## SC and I4 comparison
# load IPA pathway comparison z-score
PW_Z <- read_excel("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/I4/SCXI4_IPA.xls",skip=1)
PW_Z <- column_to_rownames(PW_Z,"Canonical Pathways")
# load IPA pathway comparison adj.p
PW_P <- read_excel("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/I4/SCXI4_IPA_adjP.xls",skip=1)
PW_P <- column_to_rownames(PW_P,"Canonical Pathways")

# drop the insignificant pathways
PW_P <- PW_P %>% filter(SC > 1.3 & I4 > 1.3)
PW_Z <- PW_Z[rownames(PW_P),]

PW_Z <- as.matrix(PW_Z)
PW_P <- as.matrix(PW_P)

# replace Greek symbols
library(stringi)
Greek <- c("α","β","γ","κ","θ")
English <- c("a","b","r","k","th")
rownames(PW_Z) <- stri_replace_all_regex(rownames(PW_Z),
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)
rownames(PW_P) <- stri_replace_all_regex(rownames(PW_P),
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)

# reorder the p value matrix by the z-score matrix
PW_P <- PW_P[rownames(PW_Z),colnames(PW_Z)]

# draw heatmap
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(file="IPA_heatmap_p.pdf", width=7, height = 14)
ht <- Heatmap(PW_Z,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              cluster_columns = F,
              column_dend_height = unit(18, "mm"),
              heatmap_legend_param = list(title="Activation\nz-score",
                                          at=c(-2,-1,0,1,2),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(PW_Z)),
              row_names_gp = gpar(fontsize = 9.5))
ht<-draw(ht)
dev.off()


## SC and Twin comparison
# load IPA pathway comparison z-score
PW_Z <- read_excel("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Twins/comparison/SCXTwin_LD.xls",skip=1)
PW_Z <- column_to_rownames(PW_Z,"Canonical Pathways")
PW_Z <- dplyr::mutate_all(PW_Z, function(x) as.numeric(as.character(x))) # convert entire datafrome to numberic
PW_Z[is.na(PW_Z)] <- 0

# load IPA pathway comparison adj.p
PW_P <- read_excel("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Twins/comparison/SCXTwin_LD_adjP.xls",skip=1)
PW_P <- column_to_rownames(PW_P,"Canonical Pathways")

# drop the insignificant pathways
PW_P <- PW_P %>% filter(SC_LD > 1.3 & Twin_LD > 1.3)
PW_Z <- PW_Z[rownames(PW_P),]

PW_Z <- as.matrix(PW_Z)
PW_P <- as.matrix(PW_P)

# replace Greek symbols
library(stringi)
Greek <- c("α","β","γ","κ","θ")
English <- c("a","b","r","k","th")
rownames(PW_Z) <- stri_replace_all_regex(rownames(PW_Z),
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)
rownames(PW_P) <- stri_replace_all_regex(rownames(PW_P),
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)

# reorder the p value matrix by the z-score matrix
PW_P <- PW_P[rownames(PW_Z),colnames(PW_Z)]

# draw heatmap
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

pdf(file="IPA_heatmap_p.pdf", width=7, height = 8)
ht <- Heatmap(PW_Z,
              col = col_fun,
              row_dend_side = "left", column_dend_side = "top",
              column_names_side = "bottom",
              cluster_columns = F,
              column_dend_height = unit(18, "mm"),
              heatmap_legend_param = list(title="Activation\nz-score",
                                          at=c(-2,-1,0,1,2),
                                          legend_gp = gpar(fontsize = 20)),
              row_names_max_width = max_text_width(rownames(PW_Z)),
              row_names_gp = gpar(fontsize = 9.5))
ht<-draw(ht)
dev.off()