read_maf <- function(session,maf_file_path) {
  maf_progress <- shiny::Progress$new(session,max=100)
  maf_progress$set(value = 0, message = "Reading MAF file...")
  maf <- read.maf(maf_file_path)

  maf_progress$set(value = 100, message = "Upload complete.")
  maf_progress$close()
  
  return(maf)  
}

make_gene_ribbon_plot <- function(session, input, plot_values) {
  
  mycoloropt <- NULL
  if (input$gene_comut_customribboncolor) {
    mycoloropt <- input$gene_comut_ribbon_color
  }
  
  plotopts <- list(maf_obj=plot_values$maf_obj, 
                   ribbon_color=mycoloropt, 
                   topN=isolate(input$gene_comut_topN),
                   pval_low=isolate(input$gene_comut_pvalLow), 
                   pval_high=isolate(input$gene_comut_pvalHigh),
                   plot_file=isolate(plot_values$gene_interaction_pdf),
                   plot_frac_mut_axis=TRUE,
                   scale_ribbon_to_fracmut=TRUE)
  
  ribbon_progress <- shiny::Progress$new(session,max=100)
  ribbon_progress$set(value = 5, message = "Making co-mutated gene ribbon plot...")
  
  updateProgress <- function(value = NULL, detail = NULL) {
    if (is.null(value)) {
      value <- ribbon_progress$getValue()
      value <- value + (ribbon_progress$getMax() - value) / 5
    }
    ribbon_progress$set(value = value, message = detail)
  }
  
  make_single_ribbon_plot(plotopts$maf_obj, 
                          ribbon_color=plotopts$ribbon_color, 
                          topN=plotopts$topN,
                          pval_low=plotopts$pval_low, 
                          pval_high=plotopts$pval_high,
                          plot_file=plotopts$plot_file,
                          progress_func=plotopts$progress_func,
                          plot_frac_mut_axis=plotopts$plot_frac_mut_axis,
                          scale_ribbon_to_fracmut=plotopts$scale_ribbon_to_fracmut)  
  ribbon_progress$set(value = 100, message = "Returning plot...")
  ribbon_progress$close()
  
}

TCGA_project_choices=c(
  "[TCGA-LAML] Acute Myeloid Leukemia"="TCGA-LAML",
  "[TCGA-ACC] Adrenocortical carcinoma"="TCGA-ACC",
  "[TCGA-BLCA] Bladder Urothelial Carcinoma"="TCGA-BLCA",
  "[TCGA-LGG] Brain Lower Grade Glioma"="TCGA-LGG",
  "[TCGA-BRCA] Breast invasive carcinoma"="TCGA-BRCA",
  "[TCGA-CESC] Cervical squamous cell carcinoma and endocervical adenocarcinoma"="TCGA-CESC",
  "[TCGA-CHOL] Cholangiocarcinoma"="TCGA-CHOL",
  "[TCGA-LCML] Chronic Myelogenous Leukemia"="TCGA-LCML",
  "[TCGA-COAD] Colon adenocarcinoma"="TCGA-COAD",
  # "[TCGA-CNTL] Controls"="TCGA-CNTL",
  "[TCGA-ESCA] Esophageal carcinoma"="TCGA-ESCA",
  # "[TCGA-FPPP] FFPE Pilot Phase II"="TCGA-FPPP",
  "[TCGA-GBM] Glioblastoma multiforme"="TCGA-GBM",
  "[TCGA-HNSC] Head and Neck squamous cell carcinoma"="TCGA-HNSC",
  "[TCGA-KICH] Kidney Chromophobe"="TCGA-KICH",
  "[TCGA-KIRC] Kidney renal clear cell carcinoma"="TCGA-KIRC",
  "[TCGA-KIRP] Kidney renal papillary cell carcinoma"="TCGA-KIRP",
  "[TCGA-LIHC] Liver hepatocellular carcinoma"="TCGA-LIHC",
  "[TCGA-LUAD] Lung adenocarcinoma"="TCGA-LUAD",
  "[TCGA-LUSC] Lung squamous cell carcinoma"="TCGA-LUSC",
  "[TCGA-DLBC] Lymphoid Neoplasm Diffuse Large B-cell Lymphoma"="TCGA-DLBC",
  "[TCGA-MESO] Mesothelioma"="TCGA-MESO",
  "[TCGA-MISC] Miscellaneous"="TCGA-MISC",
  "[TCGA-OV] Ovarian serous cystadenocarcinoma"="TCGA-OV",
  "[TCGA-PAAD] Pancreatic adenocarcinoma"="TCGA-PAAD",
  "[TCGA-PCPG] Pheochromocytoma and Paraganglioma"="TCGA-PCPG",
  "[TCGA-PRAD] Prostate adenocarcinoma"="TCGA-PRAD",
  "[TCGA-READ] Rectum adenocarcinoma"="TCGA-READ",
  "[TCGA-SARC] Sarcoma"="TCGA-SARC",
  "[TCGA-SKCM] Skin Cutaneous Melanoma"="TCGA-SKCM",
  "[TCGA-STAD] Stomach adenocarcinoma"="TCGA-STAD",
  "[TCGA-TGCT] Testicular Germ Cell Tumors"="TCGA-TGCT",
  "[TCGA-THYM] Thymoma"="TCGA-THYM",
  "[TCGA-THCA] Thyroid carcinoma"="TCGA-THCA",
  "[TCGA-UCS] Uterine Carcinosarcoma"="TCGA-UCS",
  "[TCGA-UCEC] Uterine Corpus Endometrial Carcinoma"="TCGA-UCEC",
  "[TCGA-UVM] Uveal Melanoma"="TCGA-UVM"
)


TCGA_pipeline_choices=c(
  "muse", 
  "varscan2", 
  "somaticsniper", 
  "mutect2"
)