# Prepare the PBMC uG reference sample file for CIBERSORTx

# df1 is converted of sparse matrix to a dense dataframe
df1 <- as.data.frame(Sample.combined@assays[["RNA"]]@counts)
# df2 has one column and row names
df2 <- Sample.combined@meta.data %>% select(predicted.celltype.l2)
# Create a named vector from df2 with row names as keys and the column values as values
mapping <- setNames(as.vector(df2[,1]), rownames(df2))
# Change the column names of df1 based on the mapping
colnames(df1) <- sapply(colnames(df1), function(x) ifelse(x %in% names(mapping), mapping[x], x))

# Save the dense dataframe to a tab-delimited .txt file without double quotations
write.table(df1, file = "PBMC_uG_counts_matrix.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(df1[,1:5], file = "PBMC_uG_counts_matrix_test.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)


# Use unstimulated PBMC only
# df1 is converted of sparse matrix to a dense dataframe
df1 <- as.data.frame(Sample.combined.unstim@assays[["RNA"]]@counts)
df1 <- tibble::rownames_to_column(df1, var = "Gene symbol")
# df2 has one column and row names
df2 <- Sample.combined.unstim@meta.data %>% select(predicted.celltype.l2)
# Create a named vector from df2 with row names as keys and the column values as values
mapping <- setNames(as.vector(df2[,1]), rownames(df2))
# Change the column names of df1 based on the mapping
colnames(df1) <- sapply(colnames(df1), function(x) ifelse(x %in% names(mapping), mapping[x], x))

# Save the dense dataframe to a tab-delimited .txt file without double quotations
write.table(df1, file = "PBMC_uG_unstim_counts_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# prepare to plot the results
CIBERSORTx_Results <- read.delim("cybrsortx_unstim/CIBERSORTx_Results.txt")
CIBERSORTx_Results <- CIBERSORTx_Results[1:27]
CIBERSORTx_Results4plot <- pivot_longer(CIBERSORTx_Results, cols = 2:27,names_to = "cell_type", values_to = "Freq")

## set the levels in order we want

CIBERSORTx_Results4plot$group <- case_when(CIBERSORTx_Results4plot$Mixture == 
                                             "h1G_1201_A" ~ "1G",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1G_1207_E" ~ "1G",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1G_1213_I" ~ "1G",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1G_0720" ~ "1G",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1G_0807" ~ "1G",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1G_0810" ~ "1G",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uG_1201_C" ~ "uG",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uG_1207_G" ~ "uG",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uG_1213_K" ~ "uG",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uG_0720" ~ "uG",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uG_0807" ~ "uG",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uG_0810" ~ "uG",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1GwQ_1201_B" ~ "1G_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1GwQ_1207_F" ~ "1G_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1GwQ_1213_J" ~ "1G_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1GwQ0720" ~ "1G_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1GwQ_0807" ~ "1G_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "h1GwQ_0810" ~ "1G_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uGwQ_1201_D" ~ "uG_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uGwQ_1207_H" ~ "uG_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uGwQ_1213_L" ~ "uG_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uG_wQ_0720" ~ "uG_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uGwQ_0807" ~ "uG_Q",
                                           CIBERSORTx_Results4plot$Mixture == 
                                             "uGwQ_0810" ~ "uG_Q")
                                           
CIBERSORTx_Results4plot$group<-factor(CIBERSORTx_Results4plot$group, 
                            levels=c("1G","uG","1G_Q","uG_Q"))

# All 4 groups that containing Q treated samples
# Calculate summary statistics
CIBERSORTx_Results4plot_summary <- CIBERSORTx_Results4plot %>%
  group_by(cell_type, group) %>%
  summarise(Freq_mean = mean(Freq),
            se = sd(Freq) / sqrt(n()))

# Create the bar plot with T-shaped error bars
bar_plot <- ggplot(CIBERSORTx_Results4plot_summary, aes(x = cell_type, y = Freq_mean, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Freq_mean - se, ymax = Freq_mean + se), position = position_dodge(0.9)) +
  ggtitle("Bulk RNAseq Cell Type Frequency Predicted by CIBERSORTx") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell types", y = "Proportion") +
  guides(fill = guide_legend(title = 'Groups'))

# Add the individual "Freq" values as dots
dot_plot <- ggplot(CIBERSORTx_Results4plot, aes(x = cell_type, y = Freq, fill = group)) +
  geom_jitter(size = 1, width = 0.2 / length(unique(CIBERSORTx_Results4plot$group))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell types", y = "Proportion") +
  guides(fill = guide_legend(title = 'Groups'))

# Combine the bar plot with error bars and the dot plot
combined_plot <- bar_plot + 
  geom_point(data = CIBERSORTx_Results4plot, aes(x = cell_type, y = Freq, fill = group),
             position = position_dodge(0.9), size = 0.5, alpha=0.6)
print(combined_plot)
ggsave("Cell_type_frequency_unstimu_all_CIBERSORTx.pdf",width=10,height=6)


# select 1G and uG only
sub.prop.1GuG <- subset(CIBERSORTx_Results4plot, group %in% c("1G","uG"))
sub.prop.1GuG.summary <- subset(CIBERSORTx_Results4plot_summary, group %in% c("1G","uG"))
# Create the bar plot with T-shaped error bars
bar_plot <- ggplot(sub.prop.1GuG.summary, aes(x = cell_type, y = Freq_mean, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Freq_mean - se, ymax = Freq_mean + se), position = position_dodge(0.9)) +
  ggtitle("Bulk RNAseq Cell Type Frequency Predicted by CIBERSORTx") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell types", y = "Proportion") +
  guides(fill = guide_legend(title = 'Groups'))

# Add the individual "Freq" values as dots
dot_plot <- ggplot(sub.prop.1GuG, aes(x = cell_type, y = Freq, fill = group)) +
  geom_jitter(size = 1, width = 0.2 / length(unique(sub.prop.1GuG$group))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cell types", y = "Proportion") +
  guides(fill = guide_legend(title = 'Groups'))

# Combine the bar plot with error bars and the dot plot
combined_plot <- bar_plot + 
  geom_point(data = sub.prop.1GuG, aes(x = cell_type, y = Freq, fill = group),
             position = position_dodge(0.9), size = 0.5, alpha=0.6)
print(combined_plot)
ggsave("Cell_type_frequency_unstimu_1GuG_CIBERSORTx.pdf",width=10,height=6)


# cell type ratio change: log2FC comparision unstim; remove average celltype ratio <1%
sub.prop.1GuG.1<-subset(sub.prop.1GuG.summary, group=="1G" & Freq_mean>0.01)
sub.prop.1GuG.u<-subset(sub.prop.1GuG.summary, group=="uG" & Freq_mean>0.01)
sub.prop.1GuG.comm<-intersect(sub.prop.1GuG.1$cell_type,sub.prop.1GuG.u$cell_type)
sub.prop.1GuG.1 <- subset(sub.prop.1GuG.1, cell_type %in% sub.prop.1GuG.comm)
sub.prop.1GuG.u <- subset(sub.prop.1GuG.u, cell_type %in% sub.prop.1GuG.comm)
sub.prop.1GvsuG <- data.frame(Cell_type=sub.prop.1GuG.1$cell_type,
                         group="uG vs 1G unstimulated",
                         Fold_change=sub.prop.1GuG.u$Freq_mean/sub.prop.1GuG.1$Freq_mean)
sub.prop.1GvsuG$log2FC<-log2(sub.prop.1GvsuG$Fold_change)

# add t-test p-value
sub.prop.1GuG.1.full <- subset(sub.prop.1GuG, group=="1G")
sub.prop.1GuG.u.full <- subset(sub.prop.1GuG, group=="uG")

t.test.num <- numeric()
for (c in sub.prop.1GuG.comm){
  sub.prop.1GuG.1.temp <- sub.prop.1GuG.1.full %>% filter(cell_type == c) %>% select(Freq)
  sub.prop.1GuG.u.temp <- sub.prop.1GuG.u.full %>% filter(cell_type == c) %>% select(Freq)
  t.test.temp <- t.test(sub.prop.1GuG.1.temp, sub.prop.1GuG.u.temp)
  t.test.num <- c(t.test.num,t.test.temp$p.value)
}
sub.prop.1GvsuG$p_value <- t.test.num
# bar plot
ggplot(sub.prop.1GvsuG, aes(Cell_type,log2FC))+
  geom_bar(position="dodge", stat="identity",fill="#999999") +
  ggtitle("uG vs 1G unstimulated") +
  geom_text(aes(label=round(log2FC,2)),size=2.5,
            position=position_dodge(0.9), vjust=-0.5) +
  geom_text(aes(label = ifelse(Cell_type == "CD14.Mono", "*", "")),
            position = position_dodge(0.9), vjust = 1.5,size=5) +
  theme_bw() +
  theme(plot.title = element_text(hjust =0.5),
        axis.text.x = element_text(angle=45, hjust=1))
ggsave("CellType_ratio_uG_vs_1G_unstimulated_CIBERSORTx.pdf",width=6,height=4)