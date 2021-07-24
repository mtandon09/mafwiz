
get_tcga_data <- function(tcga_dataset="ACC",save_folder=file.path("data"), variant_caller="mutect2") {
  
  require(TCGAbiolinks)
  message("Getting TCGA MAF...")
  tcga_dataset <- gsub("TCGA-","",tcga_dataset)
  save_folder_base=file.path(save_folder,paste0("TCGA-",tcga_dataset))
  save_folder=file.path(save_folder_base,variant_caller)
  tcga_maf_file=file.path(save_folder,paste0("TCGA-",tcga_dataset,".",variant_caller,".maf"))
  
  if (!file.exists(tcga_maf_file)) {
    if(!dir.exists(save_folder)) {dir.create(save_folder, recursive = T)}
    tcga_maf <- GDCquery_Maf(tcga_dataset, 
                             pipelines = variant_caller, 
                             directory = save_folder)
    tcga_maf$Tumor_Sample_Barcode_original <- tcga_maf$Tumor_Sample_Barcode
    ### The codes in the MAF are specimen IDs, but clinical data is provided by patient IDs
    ## The correct solution would be get the biospecimen data and clinical data and link it together to get patient IDs
    ## But here I'm just chopping the specimen part of the code off
    tcga_maf$Tumor_Sample_Barcode <-unlist(lapply(strsplit(tcga_maf$Tumor_Sample_Barcode, "-"), function(x) {paste0(x[1:3], collapse="-")}))
    tcga_maf$caller <- variant_caller
    write.table(tcga_maf, file=tcga_maf_file, quote=F, sep="\t", row.names = F, col.names = T)
  }
  
  
  message("Adding clinical data")
  tcga_clinical_file=file.path(save_folder_base,paste0("TCGA-",tcga_dataset,".clinical.txt"))
  if (! file.exists(tcga_clinical_file)) {
    if (!dir.exists(dirname(tcga_clinical_file))) {dir.create(dirname(tcga_clinical_file), recursive = T)}
    tcga_clinical <- GDCquery_clinic(project = paste0("TCGA-",tcga_dataset), type = "clinical")
    tcga_clinical$Tumor_Sample_Barcode <- tcga_clinical$bcr_patient_barcode
    write.table(tcga_clinical, file=tcga_clinical_file, quote=T, sep="\t", row.names = F, col.names = T)
  }
  
  tcga_clin_data <- read.table(tcga_clinical_file, sep="\t",header = T,stringsAsFactors = F)
  tcga_maf <- read.maf(tcga_maf_file, clinicalData = tcga_clin_data)
  
  # return(list(mafObj=tcga_maf, clindat=tcga_clin_data))
  return(tcga_maf)
  
}


make_tcga_clinical_annotation <- function(tcga_maf_obj, oncomat_to_match=NULL) {
  require(maftools)
  require(RColorBrewer)
  require(ComplexHeatmap)
  require(circlize)
  tcga_clin_data <- tcga_maf_obj@clinical.data
  tcga_pheno_columns <- c("Tumor_Sample_Barcode","ajcc_pathologic_stage","age_at_diagnosis","gender","race","vital_status","tissue_or_organ_of_origin")
  matched_order=1:nrow(tcga_clin_data)
  if (!is.null(oncomat_to_match)) {
    matched_order=match(colnames(oncomat_to_match), tcga_clin_data$Tumor_Sample_Barcode, nomatch=0)
  } 
  tcga_anno_data <- tcga_clin_data[matched_order,..tcga_pheno_columns]
  tcga_dataset <- paste0(unique(tcga_clin_data$disease), collapse=",")
  tcga_anno_data$Dataset <- tcga_dataset
  
  anno_data <- tcga_anno_data
  
  stages=sort(unique(anno_data$ajcc_pathologic_stage))
  stage_colors <- setNames(brewer.pal(n = length(stages), name = "Reds"), stages)
  
  anno_data$age_at_diagnosis <- as.numeric(as.character(anno_data$age_at_diagnosis))
  age_range=round(range(anno_data$age_at_diagnosis, na.rm = T),-1)
  age_color_length=10
  age_breaks=round(seq(age_range[1], age_range[2], length.out=age_color_length),0)
  age_color_vals=colorRampPalette(c("lightblue1","royalblue1","navy"))(age_color_length)
  age_colors=colorRamp2(age_breaks, age_color_vals)
  
  gender_colors=c(female="hotpink", male="cornflowerblue")
  
  races=sort(unique(anno_data$race))
  race_colors <- setNames(rev(brewer.pal(n = length(races), name = "Set1")), races)
  
  statuses=sort(unique(anno_data$vital_status))
  vitstat_colors <- c(Alive="darkgreen",Dead="darkred")
  
  tissues=sort(unique(anno_data$tissue_or_organ_of_origin))
  tissue_colors <- setNames(brewer.pal(n = length(tissues), name = "Dark2"), tissues)
  
  # dataset_colors <- setNames(c("mediumorchid1","darkolivegreen1"),
  dataset_colors <- setNames(c("grey30","darkolivegreen1"),
                             c(tcga_dataset, "Other"))
  
  anno_colors <- setNames(list(stage_colors, age_colors, gender_colors, race_colors, vitstat_colors, tissue_colors, dataset_colors),
                          setdiff(colnames(anno_data),"Tumor_Sample_Barcode"))
  
  
  mycols <- which(!colnames(anno_data) %in% c("Tumor_Sample_Barcode"))
  anno_df <- anno_data[,..mycols]
  myanno <- HeatmapAnnotation(df=anno_df,col = anno_colors)
  
  return(list(colorList=anno_colors, annodata=anno_data, HManno=myanno))
  
}

tcga_clinical_colors <- function(tcga_clin_data,
                                 preset_columns=c("Tumor_Sample_Barcode","ajcc_pathologic_stage","age_at_diagnosis","gender","race","vital_status","tissue_or_organ_of_origin")
) {
  require(maftools)
  require(RColorBrewer)
  require(ComplexHeatmap)
  require(circlize)
  
  found_columns <- intersect(colnames(tcga_clin_data), preset_columns)
  if ( ! "Tumor_Sample_Barcode" %in% found_columns ) {
    stop("Clinical data must contain a 'Tumor_Sample_Barcode' column.")
  }
  found_columns <- setdiff(found_columns, c("Tumor_Sample_Barcode"))
  if (length(found_columns) < 1) {
    stop(paste0("Clinical data must contain at least one of the columns: ", paste0(preset_columns, collapse = ", ")))
  }
  
  anno_data <- tcga_clin_data[,..found_columns]
  anno_data$Dataset <- paste0(unique(tcga_clin_data$disease), collapse=",")
  
  # browser()
  color_list <- lapply(found_columns, function(featurename) {
    
    return_val=NULL
    switch(featurename,
           "ajcc_pathologic_stage"={
             stages=sort(unique(anno_data$ajcc_pathologic_stage))
             return_val <- setNames(brewer.pal(n = length(stages), name = "Reds"), stages)
             return_val <- return_val[!is.na(names(return_val))]
           },
           "age_at_diagnosis"={
             anno_data$age_at_diagnosis <- as.numeric(as.character(anno_data$age_at_diagnosis))
             age_range=round(range(anno_data$age_at_diagnosis, na.rm = T),-1)
             age_color_length=10
             age_breaks=round(seq(age_range[1], age_range[2], length.out=age_color_length),0)
             age_color_vals=colorRampPalette(c("lightblue1","royalblue1","navy"))(age_color_length)
             return_val <- colorRamp2(age_breaks, age_color_vals)
           },
           "gender"={
             return_val <- c(female="hotpink", male="cornflowerblue")
           },
           "race"={
             races=sort(unique(anno_data$race))
             return_val <- setNames(rev(brewer.pal(n = length(races), name = "Set1")), races)
             return_val <- return_val[!is.na(names(return_val))]
           },
           "vital_status"={
             # statuses=sort(unique(anno_data$vital_status))
             return_val <- c(Alive="darkgreen",Dead="darkred")
           },
           "tissue_or_organ_of_origin"={
             tissues=sort(unique(anno_data$tissue_or_organ_of_origin))
             return_val <- setNames(brewer.pal(n = length(tissues), name = "Dark2"), tissues)
             return_val <- return_val[!is.na(names(return_val))]
           },{
             warning(paste0("'",featurename,"' not found in columns of clinical data provided."))
           }
    )
    return(return_val)
  })
  names(color_list) <- found_columns
  return(color_list)
}



get_tcga_rna_data <- function(curr_dataset, 
                              datatype="HTSeq - FPKM-UQ",
                              download_dir = "GDCdata", 
                              summarizedExperiment=T,
                              legacy=F) {
  # query <- GDCquery(project = paste0("TCGA-",gsub("TCGA-","",curr_dataset)),
  query <- GDCquery(project = curr_dataset,
                    # data.category = "Gene expression",
                    data.category = ifelse(legacy,"Gene expression","Transcriptome Profiling"),
                    workflow.type = datatype,
                    file.type  = ifelse(legacy,"results",""),
                    data.type = "Gene Expression Quantification",
                    legacy = legacy)
  # legacy = T)
  if (!dir.exists(download_dir)){dir.create(download_dir, recursive = T)}
  GDCdownload(query, method = "api", directory = download_dir)
  data <- GDCprepare(query, directory = download_dir,summarizedExperiment=summarizedExperiment)
  
  if (summarizedExperiment==F) {
    # browser()
    # query <- GDCquery(project = curr_dataset,
    # query <- GDCquery(project = c("TARGET-ALL-P1","TARGET-ALL-P2"),
    #                   data.category = "Clinical",
    #                   # data.type = "Clinical Supplement",
    #                   data.format = "XLSX",
    #                   legacy = legacy)
    # # legacy = T)
    # if (!dir.exists(download_dir)){dir.create(download_dir, recursive = T)}
    # GDCdownload(query, method = "api", directory = download_dir)
    # clinical_data <- GDCprepare(query, directory = download_dir)
    
  }
  
  
  return(data)
}

