# pepFunk

## A peptide-centric functional workflow for gut microbiome metaproteomic data

This is the source code for the [pepFunk R shiny app](https://shiny.imetalab.ca/pepFunk). 


## How to run the app

The easiest way to run the app is by using [RStudio](https://rstudio.com).
After installing the required R packages, open `app.R` in RStudio and select "Run app" at the top right corner of the RStudio window.
The app will open in a window that will allow you to run the Shiny app locally. 


## Required R packages:
- `rhandsontable`
- `shiny`
- `shinydashboard`
- `shinyWidgets`
- `colourpicker`
- `reshape2`
- `DT`
- `tidyverse`
- `plyr` 
- `DESeq2` 
- `GSVA` 
- `limma` 
- `ggdendro` 
- `plotly` 
- `dendextend` 
- `LaCroixColoR`(www.github.com/johannesbjork/LaCroixColoR)
- `shinycssloaders` 


## Tested R and package versions:

This R shiny application was tested on the following R version using the following package versions:

```
R version 3.6.2 (2019-12-12)
Platform: x86_64-apple-darwin18.7.0 (64-bit)
Running under: macOS Catalina 10.15.2

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /usr/local/Cellar/openblas/0.3.7/lib/libopenblasp-r0.3.7.dylib

locale:
[1] en_CA.UTF-8/en_CA.UTF-8/en_CA.UTF-8/C/en_CA.UTF-8/en_CA.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] shinycssloaders_0.2.0       LaCroixColoR_0.1.0          dendextend_1.13.2          
 [4] plotly_4.9.0                ggdendro_0.1-20             limma_3.42.0               
 [7] GSVA_1.32.0                 DESeq2_1.24.0               SummarizedExperiment_1.14.1
[10] DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.55.0         
[13] Biobase_2.44.0              GenomicRanges_1.36.1        GenomeInfoDb_1.20.0        
[16] IRanges_2.18.3              S4Vectors_0.22.1            BiocGenerics_0.32.0        
[19] plyr_1.8.5                  forcats_0.4.0               stringr_1.4.0              
[22] dplyr_0.8.3                 purrr_0.3.2                 readr_1.3.1                
[25] tidyr_1.0.0                 tibble_2.1.3                ggplot2_3.2.1              
[28] tidyverse_1.2.1             DT_0.10                     reshape2_1.4.3             
[31] colourpicker_1.0            shinyWidgets_0.4.9          shinydashboard_0.7.1       
[34] rhandsontable_0.3.7         shiny_1.4.0                

loaded via a namespace (and not attached):
 [1] colorspace_1.4-1       rsconnect_0.8.15       htmlTable_1.13.2       XVector_0.24.0        
 [5] base64enc_0.1-3        rstudioapi_0.10        farver_2.0.3           bit64_0.9-7           
 [9] AnnotationDbi_1.46.1   fansi_0.4.1            lubridate_1.7.4        xml2_1.2.2            
[13] splines_3.6.2          geneplotter_1.62.0     knitr_1.25             shinythemes_1.1.2     
[17] zeallot_0.1.0          Formula_1.2-3          jsonlite_1.6           packrat_0.5.0         
[21] broom_0.5.2            annotate_1.62.0        cluster_2.1.0          graph_1.62.0          
[25] compiler_3.6.2         httr_1.4.1             backports_1.1.5        assertthat_0.2.1      
[29] Matrix_1.2-18          fastmap_1.0.1          lazyeval_0.2.2         cli_2.0.1             
[33] later_1.0.0            acepack_1.4.1          htmltools_0.4.0        tools_3.6.2           
[37] gtable_0.3.0           glue_1.3.1             GenomeInfoDbData_1.2.2 Rcpp_1.0.3            
[41] cellranger_1.1.0       vctrs_0.2.1            nlme_3.1-142           xfun_0.10             
[45] rvest_0.3.4            mime_0.7               miniUI_0.1.1.1         lifecycle_0.1.0       
[49] XML_3.98-1.20          MASS_7.3-51.4          zlibbioc_1.30.0        scales_1.1.0          
[53] hms_0.5.1              promises_1.1.0         RColorBrewer_1.1-2     yaml_2.2.0            
[57] memoise_1.1.0          gridExtra_2.3          rpart_4.1-15           latticeExtra_0.6-28   
[61] stringi_1.4.3          RSQLite_2.1.2          genefilter_1.66.0      checkmate_1.9.4       
[65] rlang_0.4.2            pkgconfig_2.0.3        bitops_1.0-6           lattice_0.20-38       
[69] htmlwidgets_1.5.1      bit_1.1-14             tidyselect_0.2.5       GSEABase_1.46.0       
[73] magrittr_1.5           R6_2.4.1               generics_0.0.2         Hmisc_4.2-0           
[77] DBI_1.1.0              pillar_1.4.3           haven_2.1.1            foreign_0.8-72        
[81] withr_2.1.2            survival_3.1-8         RCurl_1.95-4.12        nnet_7.3-12           
[85] modelr_0.1.5           crayon_1.3.4           viridis_0.5.1          locfit_1.5-9.1        
[89] grid_3.6.2             readxl_1.3.1           data.table_1.12.8      blob_1.2.0            
[93] digest_0.6.23          xtable_1.8-4           httpuv_1.5.2           munsell_0.5.0         
[97] viridisLite_0.3.0 
```
