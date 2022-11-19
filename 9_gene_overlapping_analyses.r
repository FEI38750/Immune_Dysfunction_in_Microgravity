### --- Gene overlapping analyses --- ###
## SC X Bulk RNA-seq
library(dplyr)
# read SC result
SC <- read.csv("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Markers/Unstimulated/markers_uGvs1G_375_unstim.csv")
colnames(SC)[1]<-"gene_name"
SC <- SC %>% dplyr::select(gene_name,avg_log2FC,p_val,p_val_adj) %>% 
  filter(p_val_adj < 0.05)
colnames(SC)[-1] <- paste0(colnames(SC)[-1],"_SC")
# SC_up
SC_up <- SC %>% filter(avg_log2FC_SC > 0)
# SC_down
SC_down <- SC %>% filter(avg_log2FC_SC < 0)

# read bulk result
bulk<-read.csv("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Bulk_RNAseq/Host_DEG_n3/host_counts_DEG.csv")
bulk <- bulk %>% dplyr::select(gene_name,log2FoldChange.uG_vs_oneG,pvalue.uG_vs_oneG,padj.uG_vs_oneG,
                               chromosome_name,start_position,end_position,strand,gene_biotype,description,gene_length) %>%
  filter(pvalue.uG_vs_oneG < 0.05)
bulk <- bulk %>% filter(gene_name != "-")
colnames(bulk)[2:4] <- paste0(colnames(bulk)[2:4],"_bulk")
# bulk_up
bulk_up <- bulk %>% filter(log2FoldChange.uG_vs_oneG_bulk > 0)
# bulk_down
bulk_down <- bulk %>% filter(log2FoldChange.uG_vs_oneG_bulk < 0)

# overlapping by inner_join
SCXbulk_up <- inner_join(SC_up,bulk_up)
SCXbulk_down <- inner_join(SC_down,bulk_down)

# total overlapped gene
SCXbulk <- rbind(SCXbulk_up,SCXbulk_down)

# save overlapping results
write.csv(SCXbulk,"SCXbulk_Overlapped_genes.csv")

### --- Volcano plot --- ###
library(EnhancedVolcano)
# read the count matrix output from Bulk RNA-seq analysis 
bulk_full <-read.csv("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/Bulk_RNAseq/Host_DEG_n3/host_counts_DEG.csv")
bulk_full <- bulk_full %>% dplyr::select(hybrid_name,log2FoldChange.uG_vs_oneG,pvalue.uG_vs_oneG,padj.uG_vs_oneG)

# create custom key-value for colors
keyvals <- case_when(bulk_full$hybrid_name %in% SCXbulk_up$gene_name ~ "red4",
                     bulk_full$hybrid_name %in% SCXbulk_down$gene_name ~ "blue4",
                     bulk_full$log2FoldChange.uG_vs_oneG > 0 & bulk_full$pvalue.uG_vs_oneG <=0.05  ~ "red2",
                     bulk_full$log2FoldChange.uG_vs_oneG < 0 & bulk_full$pvalue.uG_vs_oneG <=0.05 ~ "royalblue",
                     bulk_full$pvalue.uG_vs_oneG > 0.05 ~ "grey30")

keyvals[is.na(keyvals)] <- "grey30"
names(keyvals)[keyvals == "red2"] <- "Genes UP"
names(keyvals)[keyvals == "royalblue"] <- "Genes DOWN"
names(keyvals)[keyvals == "red4"] <- "Overlapped Genes UP"
names(keyvals)[keyvals == "blue4"] <- "Overlapped Genes DOWN"
names(keyvals)[keyvals == 'grey30'] <- 'NS'

labcol <- ifelse(bulk_full$log2FoldChange.uG_vs_oneG > 0, "red4",
                 "blue4")
names(labcol) <- bulk_full$hybrid_name
labcol <- labcol[names(labcol) %in% SCXbulk$gene_name]

EnhancedVolcano(bulk_full,
                #lab = NA,
                lab = bulk_full$hybrid_name,
                selectLab = SCXbulk$gene_name,
                x = 'log2FoldChange.uG_vs_oneG',
                y = 'pvalue.uG_vs_oneG',
                FCcutoff = 0,
                pCutoff = 0.05,
                labSize = 2,
                labCol = labcol,
                colCustom = keyvals,
                boxedLabels = F,
                drawConnectors = T, arrowheads=F, min.segment.length=0.1,
                widthConnectors = 0.2,
                max.overlaps = 100,
                legendPosition = "right",
                title = "uG vs 1G",
                subtitle = bquote(italic("unstimulated PBMC Bulk RNAseq")))
ggsave("uGvs1G_unstimulated_Volcano_col.pdf",width=9.5,height=8)


## SC X GLDS-420:
# read bulk result
bulk<-read.csv("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/GLDS-420/host_counts_DEG.csv")
bulk <- bulk %>% dplyr::select(gene_name,log2FoldChange.Flight_vs_Control,pvalue.Flight_vs_Control,padj.Flight_vs_Control,
                               chromosome_name,start_position,end_position,strand,gene_biotype,description,gene_length) %>%
  filter(pvalue.Flight_vs_Control < 0.05)
colnames(bulk)[2:4] <- paste0(colnames(bulk)[2:4],"_bulk")

# convert mouse gene to human orthologous
# prepare a mouse to human symbol list table - genesV2
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = bulk$gene_name , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
# start convert
for (s in 1:length(bulk$gene_name)){
  if(bulk$gene_name[s] %in% genesV2[,1]){
    bulk$human_symbol[s] <- genesV2 %>% 
      filter(MGI.symbol == bulk$gene_name[s]) %>% 
      dplyr::select(HGNC.symbol)
  } else {bulk$human_symbol[s] <- "-"}
}
bulk <- bulk %>% filter(human_symbol != "-") # drop the absent symbols 
# for multiplex names, just use the uppercase of mouse symbol
for (s in 1:length(bulk$human_symbol)){
  if(length(bulk$human_symbol[[s]]) > 1){
    bulk$human_symbol[s] <- toupper(bulk$gene_name[s])
  }
}
bulk$human_symbol<-as.character(bulk$human_symbol)
bulk_mouse <- bulk # save the table
write.csv(bulk_mouse,"GLDS-420_orthologous.csv")
# process duplicated symbols
#bulk_mouse$human_symbol[duplicated(bulk_mouse$human_symbol)]
bulk <- bulk %>% group_by(human_symbol) %>% 
  summarise(log2FoldChange.Flight_vs_Control_bulk=mean(log2FoldChange.Flight_vs_Control_bulk),
            pvalue.Flight_vs_Control_bulk=mean(pvalue.Flight_vs_Control_bulk),
            padj.Flight_vs_Control_bulk=mean(padj.Flight_vs_Control_bulk))

names(bulk)[1] <- "gene_name"

# bulk_up
bulk_up <- bulk %>% filter(log2FoldChange.Flight_vs_Control_bulk > 0)
# bulk_down
bulk_down <- bulk %>% filter(log2FoldChange.Flight_vs_Control_bulk < 0)

# overlapping by inner_join
SCXbulk_up <- inner_join(SC_up,bulk_up)
SCXbulk_down <- inner_join(SC_down,bulk_down)

# total overlapped gene
SCXbulk <- rbind(SCXbulk_up,SCXbulk_down)

# save overlapping results
write.csv(SCXbulk,"SCXGL420_Overlapped_genes.csv")

bulk_full <-read.csv("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/GLDS-420/host_counts_DEG.csv")
bulk_full <- bulk_full %>% dplyr::select(hybrid_name,log2FoldChange.Flight_vs_Control,pvalue.Flight_vs_Control,padj.Flight_vs_Control)

# create another overlapping by using mouse symbol (including duplicated names)
  # bulk_up
  bulk_up_mouse <- bulk_mouse %>% filter(log2FoldChange.Flight_vs_Control_bulk > 0)
  # bulk_down
  bulk_down_mouse <- bulk_mouse %>% filter(log2FoldChange.Flight_vs_Control_bulk < 0)
  # overlapping
  SCXbulk_up_mouse <- inner_join(SC_up,bulk_up_mouse, by = c("gene_name" = "human_symbol"))
  SCXbulk_down_mouse <- inner_join(SC_down,bulk_down_mouse, by = c("gene_name" = "human_symbol"))

# total overlapped gene
SCXbulk_mouse <- rbind(SCXbulk_up_mouse,SCXbulk_down_mouse)
colnames(SCXbulk_mouse)[colnames(SCXbulk_mouse)=="gene_name.y"]<-"gene_name_mouse"

# create custom key-value for colors
keyvals <- case_when(bulk_full$hybrid_name %in% SCXbulk_up_mouse$gene_name.y ~ "red4",
                     bulk_full$hybrid_name %in% SCXbulk_down_mouse$gene_name.y ~ "blue4",
                     bulk_full$log2FoldChange.Flight_vs_Control > 0 & bulk_full$pvalue.Flight_vs_Control <=0.05  ~ "red2",
                     bulk_full$log2FoldChange.Flight_vs_Control < 0 & bulk_full$pvalue.Flight_vs_Control <=0.05 ~ "royalblue",
                     bulk_full$pvalue.Flight_vs_Control > 0.05 ~ "grey30")

keyvals[is.na(keyvals)] <- "grey30"
names(keyvals)[keyvals == "red2"] <- "Genes UP"
names(keyvals)[keyvals == "royalblue"] <- "Genes DOWN"
names(keyvals)[keyvals == "red4"] <- "Overlapped Genes UP"
names(keyvals)[keyvals == "blue4"] <- "Overlapped Genes DOWN"
names(keyvals)[keyvals == 'grey30'] <- 'NS'

labcol <- ifelse(bulk_full$log2FoldChange.Flight_vs_Control > 0, "red4",
                 "blue4")
names(labcol) <- bulk_full$hybrid_name
labcol <- labcol[names(labcol) %in% SCXbulk_mouse$gene_name_mouse]

EnhancedVolcano(bulk_full,
                #lab = NA,
                lab = bulk_full$hybrid_name,
                selectLab = SCXbulk_mouse$gene_name_mouse,
                x = 'log2FoldChange.Flight_vs_Control',
                y = 'pvalue.Flight_vs_Control',
                xlim = c(-5,4),
                ylim = c(0,19),
                FCcutoff = 0,
                pCutoff = 0.05,
                labSize = 2,
                labCol = labcol,
                colCustom = keyvals,
                boxedLabels = F,
                drawConnectors = T, arrowheads=F, min.segment.length=0.1,
                widthConnectors = 0.2,
                max.overlaps = 100,
                legendPosition = "right",
                title = "Flight vs Ground",
                subtitle = bquote(italic("Mouse Spleen Bulk RNAseq")))
ggsave("GLDS420_FlightvsGround_Volcano_col.pdf",width=9.5,height=8)


# SC X JAXA
# read JAXA results
library(readxl)
JAXA <- read_excel("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/JAXA_Twin/TGB_050_1_2_64samples_9group_totalcount_all0removed_scalingnormalized_SEM.xlsx")
JAXA <- JAXA %>% dplyr::select("Feature ID","Pre - Normalized means","Pre - Normalized SEM_by_Excel",
                               "Flight1 - Normalized means","Flight1 - Normalized SEM_by_Excel",
                               "Flight2 - Normalized means","Flight2 - Normalized SEM_by_Excel",
                               "Flight3 - Normalized means","Flight3 - Normalized SEM_by_Excel")
colnames(JAXA) <- c("gene_name","PreFlight_Normalized_means","PreFlight_SEM",
                    "Flight5d_Normalized_means","Flight5d_SEM",
                    "Flight30d_Normalized_means","Flight30d_SEM",
                    "Flight60d_Normalized_means","Flight60d_SEM")
# flight day 5 vs pre-flight
JAXA$FLT5DvsPre_log2FC <- log2(JAXA$Flight5d_Normalized_means/JAXA$PreFlight_Normalized_means)
JAXA<-na.omit(JAXA) # remove NA
JAXA$FLT5DvsPre_log2FC <- gsub(-Inf,sort(unique(JAXA$FLT5DvsPre_log2FC))[2],JAXA$FLT5DvsPre_log2FC) # replace -Inf
JAXA$FLT5DvsPre_log2FC <- gsub(Inf,sort(unique(JAXA$FLT5DvsPre_log2FC),decreasing = T)[2],JAXA$FLT5DvsPre_log2FC) # replace Inf
JAXA$FLT5DvsPre_Pvalue <- pt(-abs(((JAXA$Flight5d_Normalized_means-JAXA$PreFlight_Normalized_means)-0)/
                               sqrt(JAXA$Flight5d_SEM**2 + JAXA$PreFlight_SEM**2)),
                             6-1)
# flight day 30 vs pre-flight
JAXA$FLT30DvsPre_log2FC <- log2(JAXA$Flight30d_Normalized_means/JAXA$PreFlight_Normalized_means)
JAXA<-na.omit(JAXA) # remove NA
JAXA$FLT30DvsPre_log2FC <- gsub(-Inf,sort(unique(JAXA$FLT30DvsPre_log2FC))[2],JAXA$FLT30DvsPre_log2FC) # replace -Inf
JAXA$FLT30DvsPre_log2FC <- gsub(Inf,sort(unique(JAXA$FLT30DvsPre_log2FC),decreasing = T)[2],JAXA$FLT30DvsPre_log2FC) # replace Inf
JAXA$FLT30DvsPre_Pvalue <- pt(-abs(((JAXA$Flight30d_Normalized_means-JAXA$PreFlight_Normalized_means)-0)/
                               sqrt(JAXA$Flight30d_SEM**2 + JAXA$PreFlight_SEM**2)),
                             6-1)
# flight day 60 vs pre-flight
JAXA$FLT60DvsPre_log2FC <- log2(JAXA$Flight60d_Normalized_means/JAXA$PreFlight_Normalized_means)
JAXA<-na.omit(JAXA) # remove NA
JAXA$FLT60DvsPre_log2FC <- gsub(-Inf,sort(unique(JAXA$FLT60DvsPre_log2FC))[2],JAXA$FLT60DvsPre_log2FC) # replace -Inf
JAXA$FLT60DvsPre_log2FC <- gsub(Inf,sort(unique(JAXA$FLT60DvsPre_log2FC),decreasing = T)[2],JAXA$FLT60DvsPre_log2FC) # replace Inf
JAXA$FLT60DvsPre_Pvalue <- pt(-abs(((JAXA$Flight60d_Normalized_means-JAXA$PreFlight_Normalized_means)-0)/
                                     sqrt(JAXA$Flight60d_SEM**2 + JAXA$PreFlight_SEM**2)),
                              6-1)
# save the results for JAXA
write.csv(JAXA,"JAXA_5d_30d_60d_vs_preFlight.csv")
# save the results for JAXA
write.csv(JAXA %>% dplyr::select("gene_name","FLT5DvsPre_log2FC","FLT5DvsPre_Pvalue",
                                 "FLT30DvsPre_log2FC","FLT30DvsPre_Pvalue"),
          "JAXA_5d_30d_vs_preFlight.csv")


# run overlapping analysis for up and down-regulated genes: day5
bulk_JAXA_5D <- JAXA %>% dplyr::select(gene_name,FLT5DvsPre_log2FC,FLT5DvsPre_Pvalue) %>%
  filter(FLT5DvsPre_Pvalue < 0.05)
# bulk_JAXA_5D_up
bulk_JAXA_5D_up <- bulk_JAXA_5D %>% filter(FLT5DvsPre_log2FC > 0)
# bulk_JAXA_5D_down
bulk_JAXA_5D_down <- bulk_JAXA_5D %>% filter(FLT5DvsPre_log2FC < 0)

# overlapping by inner_join
SCXbulk_JAXA_5D_up <- inner_join(SC_up,bulk_JAXA_5D_up)
SCXbulk_JAXA_5D_down <- inner_join(SC_down,bulk_JAXA_5D_down)

# total overlapped gene
SCXbulk_JAXA_5D <- rbind(SCXbulk_JAXA_5D_up,SCXbulk_JAXA_5D_down)

# save overlapping results
write.csv(SCXbulk_JAXA_5D,"SCXbulk_JAXA_5D_Overlapped_genes.csv")

### --- Volcano plot --- ###
library(EnhancedVolcano)

bulk_full <- JAXA %>% dplyr::select("hybrid_name"=gene_name,"log2FoldChange.uG_vs_oneG"=FLT5DvsPre_log2FC,"pvalue.uG_vs_oneG"=FLT5DvsPre_Pvalue)
bulk_full$log2FoldChange.uG_vs_oneG<-as.numeric(bulk_full$log2FoldChange.uG_vs_oneG)

SCXbulk_up<-SCXbulk_JAXA_5D_up
SCXbulk_down<-SCXbulk_JAXA_5D_down

# create custom key-value for colors
keyvals <- case_when(bulk_full$hybrid_name %in% SCXbulk_up$gene_name ~ "red4",
                     bulk_full$hybrid_name %in% SCXbulk_down$gene_name ~ "blue4",
                     bulk_full$log2FoldChange.uG_vs_oneG > 0 & bulk_full$pvalue.uG_vs_oneG <=0.05  ~ "red2",
                     bulk_full$log2FoldChange.uG_vs_oneG < 0 & bulk_full$pvalue.uG_vs_oneG <=0.05 ~ "royalblue",
                     bulk_full$pvalue.uG_vs_oneG > 0.05 ~ "grey30")

keyvals[is.na(keyvals)] <- "grey30"
names(keyvals)[keyvals == "red2"] <- "Genes UP"
names(keyvals)[keyvals == "royalblue"] <- "Genes DOWN"
names(keyvals)[keyvals == "red4"] <- "Overlapped Genes UP"
names(keyvals)[keyvals == "blue4"] <- "Overlapped Genes DOWN"
names(keyvals)[keyvals == 'grey30'] <- 'NS'

labcol <- ifelse(bulk_full$log2FoldChange.uG_vs_oneG > 0, "red4",
                 "blue4")
names(labcol) <- bulk_full$hybrid_name
labcol <- labcol[names(labcol) %in% SCXbulk_JAXA_5D$gene_name]

EnhancedVolcano(bulk_full,
                #lab = NA,
                lab = bulk_full$hybrid_name,
                selectLab = SCXbulk_JAXA_5D$gene_name,
                x = 'log2FoldChange.uG_vs_oneG',
                y = 'pvalue.uG_vs_oneG',
                ylim = c(0,4),
                FCcutoff = 0,
                pCutoff = 0.05,
                labSize = 5,
                labCol = labcol,
                colCustom = keyvals,
                boxedLabels = F,
                drawConnectors = T, arrowheads=F,
                widthConnectors = 0.4,
                max.overlaps = 100,
                legendPosition = "right",
                title = "Day 5 vs pre-flight",
                subtitle = bquote(italic("JAXA Bulk RNAseq")))
ggsave("JAXA_Day5_RNAseq_Volcano_col.pdf",width=9.5,height=8)

# run overlapping analysis for up and down-regulated genes: day30
bulk_JAXA_30D <- JAXA %>% dplyr::select(gene_name,FLT30DvsPre_log2FC,FLT30DvsPre_Pvalue) %>%
  filter(FLT30DvsPre_Pvalue < 0.05)
# bulk_JAXA_30D_up
bulk_JAXA_30D_up <- bulk_JAXA_30D %>% filter(FLT30DvsPre_log2FC > 0)
# bulk_JAXA_30D_down
bulk_JAXA_30D_down <- bulk_JAXA_30D %>% filter(FLT30DvsPre_log2FC < 0)

# overlapping by inner_join
SCXbulk_JAXA_30D_up <- inner_join(SC_up,bulk_JAXA_30D_up)
SCXbulk_JAXA_30D_down <- inner_join(SC_down,bulk_JAXA_30D_down)

# total overlapped gene
SCXbulk_JAXA_30D <- rbind(SCXbulk_JAXA_30D_up,SCXbulk_JAXA_30D_down)

# save overlapping results
write.csv(SCXbulk_JAXA_30D,"SCXbulk_JAXA_30D_Overlapped_genes.csv")

### --- Volcano plot --- ###
library(EnhancedVolcano)

bulk_full <- JAXA %>% dplyr::select("hybrid_name"=gene_name,"log2FoldChange.uG_vs_oneG"=FLT30DvsPre_log2FC,"pvalue.uG_vs_oneG"=FLT30DvsPre_Pvalue)
bulk_full$log2FoldChange.uG_vs_oneG<-as.numeric(bulk_full$log2FoldChange.uG_vs_oneG)

SCXbulk_up<-SCXbulk_JAXA_30D_up
SCXbulk_down<-SCXbulk_JAXA_30D_down

# create custom key-value for colors
keyvals <- case_when(bulk_full$hybrid_name %in% SCXbulk_up$gene_name ~ "red4",
                     bulk_full$hybrid_name %in% SCXbulk_down$gene_name ~ "blue4",
                     bulk_full$log2FoldChange.uG_vs_oneG > 0 & bulk_full$pvalue.uG_vs_oneG <=0.05  ~ "red2",
                     bulk_full$log2FoldChange.uG_vs_oneG < 0 & bulk_full$pvalue.uG_vs_oneG <=0.05 ~ "royalblue",
                     bulk_full$pvalue.uG_vs_oneG > 0.05 ~ "grey30")

keyvals[is.na(keyvals)] <- "grey30"
names(keyvals)[keyvals == "red2"] <- "Genes UP"
names(keyvals)[keyvals == "royalblue"] <- "Genes DOWN"
names(keyvals)[keyvals == "red4"] <- "Overlapped Genes UP"
names(keyvals)[keyvals == "blue4"] <- "Overlapped Genes DOWN"
names(keyvals)[keyvals == 'grey30'] <- 'NS'

labcol <- ifelse(bulk_full$log2FoldChange.uG_vs_oneG > 0, "red4",
                 "blue4")
names(labcol) <- bulk_full$hybrid_name
labcol <- labcol[names(labcol) %in% SCXbulk_JAXA_30D$gene_name]

EnhancedVolcano(bulk_full,
                #lab = NA,
                lab = bulk_full$hybrid_name,
                selectLab = SCXbulk_JAXA_30D$gene_name,
                x = 'log2FoldChange.uG_vs_oneG',
                y = 'pvalue.uG_vs_oneG',
                ylim = c(0,4),
                FCcutoff = 0,
                pCutoff = 0.05,
                labSize = 5,
                labCol = labcol,
                colCustom = keyvals,
                boxedLabels = F,
                drawConnectors = T, arrowheads=F,
                widthConnectors = 0.2,
                max.overlaps = 100,
                legendPosition = "right",
                title = "Day 30 vs pre-flight",
                subtitle = bquote(italic("JAXA Bulk RNAseq")))
ggsave("JAXA_Day30_RNAseq_Volcano_col.pdf",width=9.5,height=8)


# SC X I4 overlaps
I4 <- read.csv("/Users/feiwu/My Drive/scRNAseq_analysis/PBMC/I4/Fei_genes_uGvs1G_390_unstimulated_I4_PBMCs.logFCpoint1.csv")
# Total overlapping by inner_join
SCXI4 <- inner_join(SC,I4,by = c("gene_name" = "gene"))
write.csv(SCXI4,"SCXI4_genes4IPA.csv") # save for pathway comparison

# overlapping by inner_join
SCXI4_up <- inner_join(SC_up,I4_up)
SCXI4_down <- inner_join(SC_up,I4_down)

df<-tibble::column_to_rownames(SCXI4,"gene_name")
df<-dplyr::select(df,avg_log2FC_SC,"avg_log2FC_I4"=avg_log2FC)
df<-dplyr::mutate_all(df, function(x) as.numeric(as.character(x))) # convert entire datafrome to numberic

### heatmap plot with gene color as directions ###
library(ComplexHeatmap)
library(dplyr)
# create custom key-value for colors
overlap_color <- case_when(rownames(df) %in% SCXI4_up$gene_name ~ "red4",
                           rownames(df) %in% SCXI4_down$gene_name ~ "blue4",
                           TRUE ~ "grey30")

names(overlap_color)[overlap_color == "red4"] <- "Consistent UP"
names(overlap_color)[overlap_color == "blue4"] <- "Consistent DOWN"
names(overlap_color)[overlap_color == 'grey30'] <- 'Inconsistent'

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

pdf(file="SCXI4_Heatmap.pdf", width=4, height = 40)
Heatmap(df,
        col = col_fun,
        cluster_columns = F,
        row_names_gp = gpar(col=overlap_color),
        column_names_gp = gpar(fontsize = 20),
        heatmap_legend_param = list(title="Log2FC",
                                       at=c(-1,-0.5,0,0.5,1),
                                       legend_gp = gpar(fontsize = 20)))
dev.off()

# divide to two penels for up and down genes
df.up <- df %>% filter(avg_log2FC_SC>0)
df.down <- df %>% filter(avg_log2FC_SC<0)

overlap_color.up <- case_when(rownames(df.up) %in% SCXI4_up$gene_name ~ "red4",
                                TRUE ~ "grey30")
names(overlap_color.up)[overlap_color.up == "red4"] <- "Consistent DOWN"
names(overlap_color.up)[overlap_color.up == 'grey30'] <- 'Inconsistent'

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

pdf(file="SCXI4_Heatmap_UP.pdf", width=4, height = 20)
Heatmap(df.up,
        col = col_fun,
        cluster_columns = F,
        row_names_gp = gpar(col=overlap_color.up),
        column_names_gp = gpar(fontsize = 20),
        heatmap_legend_param = list(title="Log2FC",
                                    at=c(-1,-0.5,0,0.5,1),
                                    legend_gp = gpar(fontsize = 20)))
dev.off()

# create custom key-value for colors
overlap_color.down <- case_when(rownames(df.down) %in% SCXI4_down$gene_name ~ "blue4",
                           TRUE ~ "grey30")
names(overlap_color.down)[overlap_color.down == "blue4"] <- "Consistent DOWN"
names(overlap_color.down)[overlap_color.down == 'grey30'] <- 'Inconsistent'

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
pdf(file="SCXI4_Heatmap_DOWN.pdf", width=4, height = 20)
Heatmap(df.down,
        col = col_fun,
        cluster_columns = F,
        row_names_gp = gpar(col=overlap_color.down),
        column_names_gp = gpar(fontsize = 20),
        heatmap_legend_param = list(title="Log2FC",
                                    at=c(-1,-0.5,0,0.5,1),
                                    legend_gp = gpar(fontsize = 20)))
dev.off()

# I4 overlaps up and down
colnames(I4)[1]<-"gene_name"
I4 <- I4 %>% dplyr::select(gene_name,avg_log2FC,p_val,p_val_adj) %>% 
  filter(p_val_adj < 0.05)
colnames(I4)[-1] <- paste0(colnames(I4)[-1],"_I4")
# I4_up
I4_up <- I4 %>% filter(avg_log2FC_I4 > 0)
# I4_down
I4_down <- I4 %>% filter(avg_log2FC_I4 < 0)

# overlapping by inner_join
SCXI4_up <- inner_join(SC_up,I4_up)
SCXI4_down <- inner_join(SC_down,I4_down)

# total overlapped gene
SCXI4_same_direction <- rbind(SCXI4_up,SCXI4_down)

# save overlapping results
write.csv(SCXI4_same_direction,"SCXI4_Overlapped_genes.csv")


### GeneOverlap ###
# gene overlap test by Fisher's Exact Test
#BiocManager::install("GeneOverlap")
library(GeneOverlap)

go.obj <- newGeneOverlap(SC[,1],bulk[,1],spec ="hg19.gene")
go.obj
go.obj <- testGeneOverlap(go.obj)
go.obj
print(go.obj)

go.obj.420 <- newGeneOverlap(SC[,1],bulk_mouse[,12],spec ="hg19.gene")
go.obj.420
go.obj.420 <- testGeneOverlap(go.obj.420)
go.obj.420
print(go.obj.420)

go.obj <- newGeneOverlap(SC[,1],SCXI4[,1],spec ="hg19.gene")
go.obj
go.obj <- testGeneOverlap(go.obj)
go.obj
print(go.obj)
