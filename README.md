# Immune_Dysfunction_in_Microgravity
## Code used in the paper: *Single Cell Analysis Identifies Conserved Features of Immune Dysfunction in Simulated Microgravity and Spaceflight*

R version 4.2.2 Patched (2022-11-10 r83330)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.5 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
 [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] corrplot_0.92               ggpubr_0.4.0                ggridges_0.5.3             
 [4] magrittr_2.0.3              SeuratWrappers_0.3.0        monocle3_1.2.9             
 [7] SingleCellExperiment_1.18.0 SummarizedExperiment_1.26.1 GenomicRanges_1.48.0       
[10] GenomeInfoDb_1.32.2         IRanges_2.30.0              S4Vectors_0.34.0           
[13] MatrixGenerics_1.8.1        matrixStats_0.62.0          Biobase_2.56.0             
[16] BiocGenerics_0.42.0         RColorBrewer_1.1-3          tidyr_1.2.0                
[19] DoubletFinder_2.0.3         dplyr_1.0.9                 tibble_3.1.8               
[22] future_1.27.0               glmGamPoi_1.8.0             sctransform_0.3.3          
[25] ggplot2_3.3.6               patchwork_1.1.1             sp_1.5-0                   
[28] SeuratObject_4.1.0          Seurat_4.1.1               

loaded via a namespace (and not attached):
  [1] readxl_1.4.1           backports_1.4.1        plyr_1.8.7            
  [4] igraph_1.3.4           lazyeval_0.2.2         splines_4.2.2         
  [7] listenv_0.8.0          scattermore_0.8        digest_0.6.29         
 [10] htmltools_0.5.3        viridis_0.6.2          fansi_1.0.3           
 [13] tensor_1.5             cluster_2.1.4          ROCR_1.0-11           
 [16] remotes_2.4.2          globals_0.16.0         R.utils_2.12.0        
 [19] spatstat.sparse_2.1-1  colorspace_2.0-3       ggrepel_0.9.1         
 [22] RCurl_1.98-1.8         jsonlite_1.8.0         lme4_1.1-30           
 [25] progressr_0.10.1       spatstat.data_2.2-0    survival_3.4-0        
 [28] zoo_1.8-10             glue_1.6.2             polyclip_1.10-0       
 [31] gtable_0.3.0           zlibbioc_1.42.0        XVector_0.36.0        
 [34] leiden_0.4.2           DelayedArray_0.22.0    car_3.1-0             
 [37] future.apply_1.9.0     abind_1.4-5            scales_1.2.0          
 [40] DBI_1.1.3              rstatix_0.7.0          spatstat.random_2.2-0 
 [43] miniUI_0.1.1.1         Rcpp_1.0.9             viridisLite_0.4.0     
 [46] xtable_1.8-4           reticulate_1.25        spatstat.core_2.4-4   
 [49] rsvd_1.0.5             htmlwidgets_1.5.4      httr_1.4.3            
 [52] ellipsis_0.3.2         ica_1.0-3              R.methodsS3_1.8.2     
 [55] pkgconfig_2.0.3        uwot_0.1.11            deldir_1.0-6          
 [58] utf8_1.2.2             tidyselect_1.1.2       rlang_1.0.4           
 [61] reshape2_1.4.4         later_1.3.0            munsell_0.5.0         
 [64] cellranger_1.1.0       tools_4.2.2            cli_3.3.0             
 [67] generics_0.1.3         broom_1.0.0            stringr_1.4.0         
 [70] fastmap_1.1.0          goftest_1.2-3          fitdistrplus_1.1-8    
 [73] purrr_0.3.4            RANN_2.6.1             pbapply_1.5-0         
 [76] nlme_3.1-161           mime_0.12              R.oo_1.25.0           
 [79] compiler_4.2.2         rstudioapi_0.13        plotly_4.10.0         
 [82] png_0.1-7              ggsignif_0.6.3         spatstat.utils_2.3-1  
 [85] stringi_1.7.8          rgeos_0.5-9            lattice_0.20-45       
 [88] Matrix_1.5-1           nloptr_2.0.3           vctrs_0.4.1           
 [91] pillar_1.8.0           lifecycle_1.0.1        BiocManager_1.30.18   
 [94] spatstat.geom_2.4-0    lmtest_0.9-40          RcppAnnoy_0.0.19      
 [97] data.table_1.14.2      cowplot_1.1.1          bitops_1.0-7          
[100] irlba_2.3.5            httpuv_1.6.5           R6_2.5.1              
[103] promises_1.2.0.1       KernSmooth_2.23-20     gridExtra_2.3         
[106] parallelly_1.32.1      codetools_0.2-18       boot_1.3-28           
[109] MASS_7.3-58            assertthat_0.2.1       withr_2.5.0           
[112] GenomeInfoDbData_1.2.8 mgcv_1.8-41            parallel_4.2.2        
[115] terra_1.6-3            grid_4.2.2             rpart_4.1.19          
[118] minqa_1.2.4            carData_3.0-5          Rtsne_0.16            
[121] shiny_1.7.2           

MTD software and its instruction is on https://github.com/FEI38750/MTD<br>
MTD is running under Conda environments with Bash, the versions of dependencies are in the .yml files in the MTD/Installation.
