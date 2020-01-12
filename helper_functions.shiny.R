read_maf <- function(session,maf_file_path) {
  maf_progress <- shiny::Progress$new(session,max=100)
  maf_progress$set(value = 0, message = "Reading MAF file...")
  maf <- read.maf(maf_file_path)

  maf_progress$set(value = 100, message = "Upload complete.")
  maf_progress$close()
  
  return(maf)  
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
  "[TCGA-CNTL] Controls"="TCGA-CNTL",
  "[TCGA-ESCA] Esophageal carcinoma"="TCGA-ESCA",
  "[TCGA-FPPP] FFPE Pilot Phase II"="TCGA-FPPP",
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