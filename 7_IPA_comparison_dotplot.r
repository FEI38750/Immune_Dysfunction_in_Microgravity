library(readxl)
library(tidyr)
library(colorspace)
#### draw dot plot for IPA comparison results
setwd("~/scRNAseq_analysis/IPA")
# read the Excel exported from IPA comparison
IPA_Z <- read_excel("unstim/unstim_IPA_selected.xls",skip=1)
IPA_P <- read_excel("unstim/unstim_IPA_selected_adjP.xls",skip=1)

# replace Greek symbols
library(stringi)
Greek <- c("α","β","γ","κ","θ")
English <- c("a","b","r","k","th")
IPA_Z$`Canonical Pathways` <- stri_replace_all_regex(IPA_Z$`Canonical Pathways`,
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)
IPA_P$`Canonical Pathways` <- stri_replace_all_regex(IPA_P$`Canonical Pathways`,
                                         pattern=Greek,
                                         replacement = English,
                                         vectorize=F)
# Order data frame rows according to vector with specific order
IPA_P<-IPA_P[match(IPA_Z$'Canonical Pathways',IPA_P$'Canonical Pathways'),]
IPA_Z <- pivot_longer(IPA_Z, cols= 2:21,names_to = "Cell_types",values_to = "Zscore")
IPA_P <- pivot_longer(IPA_P, cols= 2:21,names_to = "Cell_types",values_to = "-log(P)")
IPA_Z <- cbind(IPA_Z,`-log(P)`=IPA_P$`-log(P)`)

# trim the name of cell types
#IPA_Z$Cell_types<- gsub("_uGvs1G|_unstimulated","",IPA_Z$Cell_types)
IPA_Z$Cell_types<- gsub("_uGvs1G.*","",IPA_Z$Cell_types)
pw<-IPA_Z$'Canonical Pathways' %>% unique()
IPA_Z$Cell_types<-factor(IPA_Z$Cell_types, levels = unique(IPA_Z$Cell_types)) # order the x axis
# Replacing character values with NA in a data frame
IPA_Z[ IPA_Z == "N/A" ] <- NA
IPA_Z$Zscore <- as.numeric(IPA_Z$Zscore)

IPA_Z %>% filter(`Canonical Pathways` %in% pw) %>% filter(`-log(P)`>1.3) %>%
  ggplot(aes(x=Cell_types, y = `Canonical Pathways`, color = Zscore, size = `-log(P)`)) + 
  geom_point() +
  scale_y_discrete(labels=function(y) str_wrap(y, width=80)) +
  ylab('Canonical Pathways') +
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(size=9.5, angle=45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size=7.5))+
  theme(axis.line  = element_blank()) +
  theme(axis.ticks = element_blank()) +
  scale_color_continuous_divergingx('RdBu',rev=T,limits = c(-3,3), oob = scales::squish, name = 'z-score',
                                    na.value="transparent") +
  labs(size="-log10(adj.P)")
ggsave("Dot_uGvs1G_unstim_IPA.pdf", width = 10, height = 9,limitsize = FALSE)
