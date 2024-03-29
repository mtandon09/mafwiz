# MAF Wiz: A tool for visualizing MAF files
Maf Wiz *like math wiz*
 
## Shiny app for plotting data from MAF files
Add more docs here
 
## How to run
### Install R and RStudio
Here's an [installation guide](http://www.sthda.com/english/wiki/installing-r-and-rstudio-easy-r-programming) if you don't have these already.

### Install Required Packages
```
list.of.packages <- c("shiny","shinyjs","openxlsx","shinycssloaders","shinyWidgets","shinythemes","RColorBrewer","ggplot2","reshape2","colourpicker")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

list.of.bioc.packages <- c("maftools","ComplexHeatmap","circlize","dplyr","ggbeeswarm","TCGAbiolinks","NMF","MutationalPatterns","BSgenome.Hsapiens.UCSC.hg38","BSgenome.Hsapiens.UCSC.hg19")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
new.bioc.packages <- list.of.bioc.packages[!(list.of.bioc.packages %in% installed.packages()[,"Package"])]
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages)
```

### Run from GitHub in RStudio
```
library(shiny)
runGitHub( "mafwiz", "mtandon09", ref="master")
```

# Example usage:
Upload your own MAF file or download TCGA MAFs on-the-fly
   ![alt text](screenshots/landing.png "Select MAF")
Add sample information for annotations
   ![alt text](screenshots/sample_info.png "Add sample annotations")
Visualizations
   - Mutational Burden
    ![alt text](screenshots/burden.png "Burden Plot")
   - Oncoplot
    ![alt text](screenshots/oncoplot.png "Burden Plot")
   - Mutational Signatures (COSMIC v3)
    ![alt text](screenshots/mutational_signatures.png "Burden Plot")
   - Co-occurence of gene mutations
    ![alt text](screenshots/comut_genes.png "Burden Plot")


# To Do
- Buncha stuff
