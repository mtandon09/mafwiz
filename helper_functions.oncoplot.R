
filter_maf <- function(maf_file, flag_genes="default",save_name=NULL,no_filter=F,
                       norm_alt_max=1,t_alt_min=1,t_depth_min=20, n_callers=2,
                       gnomAD_AF_max=0.001, AF_max=0.001, ExAC_AF_max=0.001, 
                       tumor_freq_min=0.05, norm_freq_max=0.02,
                       variant_caller=NULL) {
  # browser()
  if (length(flag_genes)==0) {
    flag_genes <- c()
  } else if (flag_genes[1]=="default") {
    flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
  }
  maf_df.raw <- read.table(maf_file, sep="\t", header=T, fill = T, quote="\"", stringsAsFactors = F)
  maf_df.raw <- maf_df.raw[maf_df.raw$Hugo_Symbol != "Hugo_Symbol",]
  
  if (!"tumor_freq" %in% colnames(maf_df.raw)) {
    maf_df.raw$tumor_freq <- as.numeric(maf_df.raw$t_alt_count)/as.numeric(maf_df.raw$t_depth)
  }  
  if (!"norm_freq" %in% colnames(maf_df.raw)) {
      maf_df.raw$norm_freq <- as.numeric(maf_df.raw$n_alt_count)/as.numeric(maf_df.raw$n_depth)
  }
  
  filter_tumor_depth=rep(TRUE,nrow(maf_df.raw))
  filter_norm_alt=rep(TRUE,nrow(maf_df.raw))
  filter_tumor_alt=rep(TRUE,nrow(maf_df.raw))
  filter_genes=rep(TRUE,nrow(maf_df.raw))
  filter_pop_freq=rep(TRUE,nrow(maf_df.raw))
  
  if (!no_filter) {
    options(warn=-1)
    filter_tumor_depth=as.numeric(maf_df.raw$t_depth) > t_depth_min
    if (!sum(is.na(maf_df.raw$norm_freq)) == nrow(maf_df.raw)){
      filter_norm_alt=maf_df.raw$norm_freq < norm_freq_max
    } 
    filter_tumor_alt=maf_df.raw$tumor_freq > tumor_freq_min
    if (! is.null(t_alt_min)){
      filter_tumor_alt <- filter_tumor_alt & maf_df.raw$t_alt_count > t_alt_min
    }
    filter_genes=!maf_df.raw$Hugo_Symbol %in% flag_genes
    filter_pop_freq=(maf_df.raw$gnomAD_AF %in% c("-","") | is.na(maf_df.raw$gnomAD_AF) | as.numeric(maf_df.raw$gnomAD_AF) < min(gnomAD_AF_max,1)) & 
                    (maf_df.raw$AF %in% c("-","") | is.na(maf_df.raw$AF)  | as.numeric(maf_df.raw$AF) < min(AF_max,1)) &
                    (maf_df.raw$ExAC_AF %in% c("-","") | is.na(maf_df.raw$ExAC_AF) | as.numeric(maf_df.raw$ExAC_AF) < min(ExAC_AF_max,1))
    options(warn=0)
  }
  filter_caller=rep(TRUE,nrow(maf_df.raw))
  if (! is.null(variant_caller)) {       ### Set 'variant_caller' to NULL to skip any filtering based on caller
    maf_df.raw$set[maf_df.raw$set=="" & maf_df.raw$Hugo_Symbol=="Hugo_Symbol"] <- "set"
    maf_df.raw$set[maf_df.raw$set==""] <- "N.A."
    if (variant_caller == "consensus") {   ### Set 'variant_caller' to 'consensus' to keep variants by two or more callers
      # filter_caller <- grepl("-|Intersection", maf_df.raw$set)
      filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {length(x)>=n_callers | "Intersection" %in% x}))
    } else {                             ### Set 'variant_caller' to one of the callers (mutect, mutect2, vardict, or strelka) to get only that caller
      # filter_caller <- grepl(paste0(variant_caller,"[|-]|Intersection"), maf_df.raw$set)
      filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {any(c(variant_caller,"Intersection") %in% x)}))
    }
  }
  
  # maf_df.raw <- maf_df.raw[filter_genes &filter_norm_alt & filter_tumor_alt & filter_pop_freq & filter_caller,]
  maf_df.rawest <- maf_df.raw
  # maf_df.raw <- maf_df.rawest
  maf_df.raw <- maf_df.raw[filter_tumor_depth & filter_norm_alt & filter_tumor_alt & filter_genes & filter_pop_freq & filter_caller,]
  # browser()
  maf_df.raw <- maf_df.raw[rowSums(is.na(maf_df.raw))!=ncol(maf_df.raw),]
  # maf_df.raw <- maf_df.raw[rowSums(apply(maf_df.raw,1,is.na)]
  # maf_df.raw <- maf_df.raw[complete.cases(maf_df.raw),]
  
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    write.table(maf_df.raw, sep="\t", quote=F, file = save_name, row.names = F)
    print(paste0("Saving filtered maf to ",save_name))
    return(save_name)
  } else {
    return(maf_df.raw)
  }
}



make_oncoplot <- function(maf.filtered, cohort_freq_thresh = 0.1) {
  
  ### Read in MAF file
  # maf.filtered <- read.maf(maf_file)
  
  ### Structure info about the fraction of the cohort that has each gene mutated
  frac_mut <- data.frame(Hugo_Symbol=maf.filtered@gene.summary$Hugo_Symbol,
                         frac_mut=(maf.filtered@gene.summary$MutatedSamples/as.numeric(maf.filtered@summary$summary[3])),
                         stringsAsFactors = F)
  
  target_frac = sort(frac_mut$frac_mut, decreasing = T)[min(50,nrow(frac_mut))]*1.01
  cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
  ### Select genes based on the frequency threshold
  freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut > cohort_freq_thresh]
  if (length(freq_genes) == 0) {
    stop("No genes to plot; change the frequency threshold to include more genes.")
  }
  if (length(freq_genes) > 100) {
    target_frac = round(sort(frac_mut$frac_mut, decreasing = T)[min(50,nrow(frac_mut))],2)
    stop(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
    # return(NA)
  }
  gene_list <- list(freq_genes)
  reasons <- paste0("Cohort Freq > ",cohort_freq_thresh)
  
  ### Collect genes to plot
  genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), reason=c())
  for (i in 1:length(gene_list)) {
    if (is.na(gene_list[[i]][1])) {
      next
    }
    genes_for_oncoplot <- rbind(genes_for_oncoplot,
                                data.frame(Hugo_Symbol=gene_list[[i]],
                                           reason=reasons[i]))
  }
  genes_for_oncoplot <- cbind(genes_for_oncoplot,
                              frac=frac_mut$frac_mut[match(genes_for_oncoplot$Hugo_Symbol, frac_mut$Hugo_Symbol)])
  
  genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$reason, -genes_for_oncoplot$frac),]
  
  ### Split the oncoplot based on the reason for picking the gene
  ###   Here, we're only picked based on the frequency
  ###   But this framework is useful for plotting genes picked using various criteria
  split_idx=genes_for_oncoplot$reason
  split_colors <- rainbow(length(levels(split_idx)))
  names(split_colors) <- as.character(genes_for_oncoplot$reason[!duplicated(genes_for_oncoplot$reason)])
  split_colors <- list(Reason=split_colors)
  
  # source("scripts/helper_functions.oncoplot.R")
  ### Make matrix to plot, and order it correctly
  oncomat <- createOncoMatrix(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = F)$oncoMatrix
  oncomat <- oncomat[match(genes_for_oncoplot$Hugo_Symbol,rownames(oncomat)), ]
  onco_genes <- rownames(oncomat)
  
  oncomat.plot <- oncomat
  
  ### Set the height of the plot based on number of genes
  onco_height=NULL
  if (is.null(onco_height)) {
    onco_height=max(round(0.2*nrow(oncomat.plot),0),5)
  }
  
  ### Make the mutation type names prettier by removing the underscore
  # my_mut_col <- mutation_colors
  # names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
  oncomat.plot <- gsub("_"," ",oncomat.plot)
  
  ### Column labels get cluttered if too many samples
  show_sample_names=T
  if (ncol(oncomat.plot) > 20) {
    show_sample_names=F
  }
  # print(oncomat.plot)
  ### Make the oncoplot
  onco_base_default <- oncoPrint(oncomat.plot, alter_fun = alter_fun, col=mutation_colors, row_order=1:nrow(oncomat.plot),
                                 name="oncoplot",
                                 show_pct = F,
                                 row_split=split_idx,
                                 row_title = NULL,
                                 show_column_names = show_sample_names)#,
                                 # column_names_rot = 30,
                                 # column_gap = unit(0.0001,"npc"),
                                 # width = unit(0.75, "npc"))
  
  ### Save the oncoplot
  return(onco_base_default)
  # save_name <- paste0(out_dir,"/oncoplot.",cohort_freq_thresh,".pdf")
  # onco_width=10
  # pdf(file = save_name,height=onco_height,width=onco_width)
  # draw(onco_base_default)
  # dev.off()
}


make_burden_plot <- function(maf.filtered, plotType="Barplot") {
  num_var_data <- maf.filtered@variants.per.sample
  colnames(num_var_data) <- c("Tumor_Sample_Barcode","Variants_filtered")
  num_var_data$mut_burden <- num_var_data$Variants_filtered
  
  ## Re-jigger the factor levels so they're ordered by decreasing mutation burden (for the plotting)
  num_var_data$Tumor_Sample_Barcode <- factor(num_var_data$Tumor_Sample_Barcode,
                                              levels=num_var_data$Tumor_Sample_Barcode[order(num_var_data$mut_burden, decreasing = T)])
  
  ## Write data to a file for external plotting if desired
  # write.table(num_var_data, file = paste0(out_dir,"/mutation_burden.data.txt"), sep="\t", quote=F,row.names = F)
  
  ########################################################
  #### 5. Generate plots for mutation burden
  
  ## Pick colors
  median_mut_burdens <- num_var_data %>% summarise(median=median(mut_burden))
  
  num_var_data$xlabel <- gsub("Sample.*_","",num_var_data$Tumor_Sample_Barcode)
  num_var_data$xlabel <- factor(num_var_data$xlabel,
                                levels=num_var_data$xlabel[order(num_var_data$mut_burden, decreasing = T)])
  
  ### Mutation burden stacked with variant classification counts
  ### Works better for smaller cohorts, say < 20
  variant_type_per_sample <- as.data.frame(maf.filtered@variant.classification.summary)
  variant_type_per_sample$Tumor_Sample_Barcode <- gsub("Sample.*_","",variant_type_per_sample$Tumor_Sample_Barcode)
  var_type.melt <- reshape2::melt(variant_type_per_sample, id.vars="Tumor_Sample_Barcode",variable.name="classification",value.name="mutation_count")
  var_type.melt$mut_burden <- var_type.melt$mutation_count
  median_mut_burdens <- data.frame(median=median(var_type.melt[var_type.melt$classification== "total","mut_burden"]))
  
  plotdata <- var_type.melt[var_type.melt$classification != "total",]
  plotdata$Tumor_Sample_Barcode <- factor(as.character(plotdata$Tumor_Sample_Barcode),
                                          levels=variant_type_per_sample$Tumor_Sample_Barcode[order(variant_type_per_sample$total, decreasing = T)])
  plotdata$classification <- gsub("_"," ",plotdata$classification)
  
  class_means <- plotdata %>% group_by(classification) %>% summarise(mean=mean(mut_burden))
  plotdata$classification <- factor(as.character(plotdata$classification),
                                    levels=class_means$classification[order(class_means$mean, decreasing = F)])
  
  my_class_colors <- mutation_colors
  # names(my_class_colors) <- gsub("_", " ",names(my_class_colors))
  
  
  if (plotType=="Barplot") {
    burden_plot <- ggplot(plotdata, aes(x=Tumor_Sample_Barcode, y=mut_burden)) +
      geom_bar(aes(fill=classification), stat="identity",width=1,size=0.3, color="black") +
      scale_fill_manual(values=my_class_colors) +
      theme_linedraw(base_size = 12) +
      xlab("") + ylab("Mutations") +
      geom_hline(data = median_mut_burdens, aes(yintercept=median),linetype="dashed", color="grey60") +
      theme(
        axis.text.x = element_text(angle=30, hjust=1),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=8),
        # legend.title = element_text(size=6),
        legend.title = element_blank(),
        legend.key.height = unit(0.01,"npc"),
        legend.key.width =  unit(0.02,"npc"),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank())
  } else {
  
  
    ### Mutation Burden - Scatter/Dot plot
    ### Works better for larger cohorts
    alpha_val=1
    point_cex=2
    if (nrow(num_var_data) > 200) {
      alpha_val=0.5
    } else if (nrow(num_var_data) > 20) {
      alpha_val=0.8
    }
    burden_plot <- ggplot(num_var_data, aes(x=1, y=mut_burden)) +
      # geom_jitter(aes(x=1, y=mut_burden),inherit.aes=F,color="blue3",
      #            cex=2,size=1.5, alpha=alpha_val) +
      # geom_point(aes(x=Tumor_Sample_Barcode, y=mut_burden),inherit.aes=F,color="blue3",
      #            cex=2,size=1.5, position = position_dodge(width=0.6), alpha=alpha_val) +
      geom_beeswarm(color="blue3",cex=2,size=5,dodge.width=0.2,priority='density', alpha=alpha_val) +
      # geom_quasirandom(color="blue3",cex=2,dodge.width=1,method="tukeyDense", alpha=alpha_val) +  ### This looks nicer for fewer points/smaller cohorts than geom_beeswarm
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                   geom = "crossbar", width = 0.7, color="gray70", size = 0.2) +
      scale_y_log10()+
      theme_linedraw(base_size = 12) +
      ggtitle("Mutation Burden") +
      ylab("Mutations") + xlab("") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  return(burden_plot)
  
}



detect_maf_genome <- function(maf) {
  if (! "NCBI_Build" %in% colnames(maf@data)) {
    warning("No genome information in MAF obj.")
    return(NA)
  }
  
  my_genome = unique(maf@data$NCBI_Build)
  if (length(my_genome) > 1) {
    warning("Multiple genomes listed in MAF obj. Trying the first one")
    my_genome <- my_genome[1]
  }
  
  return_genome <- switch(my_genome,GRCh38="hg38",GRCh37="hg19",GRCm38="mm10", NA)
  
  my_chrs <- unique(maf@data$Chromosome)
  add_chr = sum(grepl("^chr", my_chrs)) < length(my_chrs)

  return(list(genome=return_genome, add_chr=add_chr))
  
}

make_mut_signature_heatmap <- function(mymaf,use_silent_mutations=F, clinVarNames = NULL, 
                                       data_dir="data", full_output=F,
                                       genome_build="hg19",
                                       progress_func=NULL) {
  
  genome_build_res <- detect_maf_genome(mymaf)
  genome_build <- genome_build_res[[1]]
  
  genome_package=paste0("BSgenome.Hsapiens.UCSC.",genome_build)
  require(genome_package, quietly = TRUE, character.only = T)
  
  # mymaf <- all_mafs[[1]]
  if (is.function(progress_func)) {
    progress_func(value=0, detail = "Computing signatures")
  }
  
  print(paste0("add chr: ",genome_build_res[[2]]))
  prefix_value=ifelse(genome_build_res[[2]],"chr","")
  # add_prefix=genome_build!="hg19"
  
  # print(date())
  print(paste0("build: ",genome_build))
  print(paste0("prefix: ",prefix_value))
  # print(add_prefix)
  print(unique(mymaf@data$Chromosome))
  # print(unique(mymaf@data$Chromosome))
  
  tnm = trinucleotideMatrix(maf = mymaf,
                            ref_genome = genome_package,
                            # add = add_prefix,
                            prefix = prefix_value,
                            useSyn = use_silent_mutations)
  mut_mat = t(tnm$nmf_matrix)
  
  # mut_mat <- data.frame(Tumor_Sample_Barcode=rownames(mut_mat), mut_mat, check.names = F)
  # mut_mat <- merge.data.frame(sample_info.exome, mut_mat, by="Tumor_Sample_Barcode")
  # write.table(mut_mat, file = paste0(figures_folder,"mutation_burden.data.txt"), sep="\t", quote=F,row.names = F)
  
  
  sp_url <- paste0(data_dir,"/cosmic/sigProfiler_exome_SBS_signatures.csv")
  cosmic_signatures = read.table(sp_url, sep = ",", header = TRUE, stringsAsFactors = F)
  cosmic_signatures$Somatic.Mutation.Type <- paste0(substr(cosmic_signatures$SubType, 1, 1),
                                                    "[",cosmic_signatures$Type, "]",
                                                    substr(cosmic_signatures$SubType, 3, 3))
  
  
  
  # Match the order of the mutation types to MutationalPatterns standard
  new_order = match(row.names(mut_mat), cosmic_signatures$Somatic.Mutation.Type)
  # Reorder cancer signatures dataframe
  cosmic_signatures = cosmic_signatures[as.vector(new_order),]
  # Add trinucletiode changes names as row.names
  row.names(cosmic_signatures) = cosmic_signatures$Somatic.Mutation.Type
  # Keep only 96 contributions of the signatures in matrix
  cosmic_signatures = as.matrix(cosmic_signatures[,grep("SBS*", colnames(cosmic_signatures))])
  
  hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
  # store signatures in new order
  cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]
  # plot(hclust_cosmic)
  
  if (is.function(progress_func)) {
    progress_func(value=50, detail = "Computing similarity to COSMIC...")
  }
  cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cosmic_signatures)
  
  # fit_res <- fit_to_signatures(mut_mat, cosmic_signatures)
  
  etio_data_file = paste0(data_dir, "/cosmic/COSMIC_signature_etiology.xlsx")
  etiology_data_raw <- read.xlsx(etio_data_file, sheet="final categories")
  etiology_data <- as.character(etiology_data_raw$CATEGORY)
  names(etiology_data) <- etiology_data_raw$signature
  etiology_colors <-  list(etiology=c("APOBEC" = "#fce116",
                                      "Defective DNA Mismatch Repair" = "#31A354",
                                      "Defective DNA Repair" = "#A1D99B",
                                      "Exonuclease Domain" = "#E5F5E0",
                                      "Exposure to Alfatoxin" = "#DE2D26",
                                      "Exposure to Aristolochic Acid" = "#FC9272",
                                      "Exposure to Haloalkanes" = "#FEE0D2",
                                      "Tobacco - Chewing" = "#6d3617",
                                      "Tobacco - Smoking" = "#a85423",
                                      "Tobacco - Smoking Associated" = "#d87841",
                                      "Prior Therapy - Alkylating Agents" = "#2171B5",
                                      "Prior Therapy - Platinum Drugs" = "#6BAED6",
                                      "Prior Therapy - Immunosuppression" = "#BDD7E7",
                                      "Prior Therapy Associated" = "#EFF3FF",
                                      "ROS Damage" = "#BCBDDC",
                                      "UV Light" = "#756BB1",
                                      "UV Light Associated" = "#a29bca",
                                      "Unknown" = "grey70",
                                      "Possible Sequencing Artifact" = "grey50")
  )

  plot_matrix <- t(cos_sim_samples_signatures)
  
  if (is.function(progress_func)) {
    progress_func(value=80, detail = "Making heatmap")
  }
  
  # browser()
  sample_info <- as.data.frame(mymaf@clinical.data)
  if (ncol(sample_info)>1) {
    sample_info <- sample_info[sample_info$Tumor_Sample_Barcode %in% colnames(plot_matrix),]
    if (is.null(clinVarNames)) {
      clinVarNames <- colnames(sample_info)
    }
    anno_data <- sample_info[, colnames(sample_info) %in% names(clinVarNames)]
    pheno_colors <- my_oncoplot_colors(anno_data)
    
  
    patient_anno <- columnAnnotation(df=anno_data, name="Patient Anno", col=pheno_colors, annotation_height = unit(1,"inches"))
    names(patient_anno) <- anno_names[names(patient_anno)]
  } else {
    mylabels <- gsub("^Sample.*?_(.*)$","\\1",sample_info$Tumor_Sample_Barcode)
    patient_anno <- columnAnnotation(Sample=anno_text(mylabels, rot = 30, gp = gpar(fontsize=10)), name="Patient Anno", annotation_height = unit(1,"inches"))
  }
  
  
  # browser()
  rowOrder=order(etiology_data)
  etiology_data <- etiology_data[rowOrder]
  signature_anno <- rowAnnotation(df=data.frame(etiology=etiology_data, row.names=rownames(plot_matrix)), 
                                  name="Signature Anno", col=etiology_colors, show_annotation_name = FALSE)
  
  # pdf(file = paste0(figures_folder,"/cosmic_v3_cosine_sim.pdf"), width=8, height=6)
  myHM <- Heatmap(plot_matrix, 
          col=colorRamp2(seq(min(plot_matrix), max(plot_matrix), length.out = 20),colorRampPalette(brewer.pal(9, "BuGn"))(20)),
          # col=colorRamp2(c(0, 0.1,0.8,1),colorRampPalette(brewer.pal(9, "YlGnBu"))(4)),
          left_annotation = signature_anno,
          # bottom_annotation = patient_anno,
          cluster_rows = F, row_order = rowOrder,
          clustering_method_rows = "median",
          clustering_method_columns = "median",
          heatmap_height = unit(6, "inches"),
          heatmap_legend_param = list(
            # at = c(-2, 0, 2),
            # labels = c("low", "zero", "high"),
            title = "Cosine Similarity",
            # legend_height = unit(4, "cm"),
            legend_direction = "horizontal"
          ),
          show_row_names=T, row_names_gp = gpar(fontsize = 5),
          show_column_names = F)
  # dev.off()
  
  if (full_output) {
    output_list <- list(plot_matrix,signature_anno, patient_anno,etiology_data, myHM)
    return(output_list)
  } else {
    return(myHM)
  }
  
}






### Cretaes matrix for oncoplot from maf file
### Adapted from maftools: https://github.com/PoisonAlien/maftools/blob/master/R/oncomatrix.R
createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){
  
  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }
  
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)
  
  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      vc = c("")
      names(vc) = 0
      
      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }
  
  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }
  
  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]
                                
                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }
                                
                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
  
  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])
  
  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)
  
  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  
  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  
  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
  
  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  
  
  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId
    
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]
    
    mdf = mdf[, -ncol(mdf)]
    
    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
    
    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
    
    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy
    
    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}

### Ordering function to make rainfall pattern
orderByGroup <- function(my_oncomat, info_vec) {
  # browser()
  class2count <- function(class_str) {
    returnval = nchar(as.character(class_str))
    if (is.na(returnval)) {
      returnval=0
    } else {
      returnval=ifelse(returnval > 2, 1, 0)
    }
    return(returnval)
  }
  oncoprint_column_order = function(count_matrix) {
    scoreCol = function(x) {
      score = 0
      for(i in 1:length(x)) {
        if(x[i]) {
          score = score + 2^(length(x)-i*1/x[i])
          score=min(score,1e16)
        }
      }
      return(score)
    }
    scores = apply(count_matrix, 2, scoreCol)
    new_order = order(scores, decreasing=TRUE)
    return(new_order)
  }
  
  oncomat_bin <- apply(my_oncomat, 1:2, class2count)
  new_names_order <- c()
  prev_max_idx=0
  for (curr_level in levels(info_vec)) {
    oncdat=oncomat_bin[,which(as.character(info_vec)==curr_level)]
    curr_new_names_order <- colnames(oncdat)[oncoprint_column_order(oncdat)]
    new_names_order <- c(new_names_order, curr_new_names_order)
  }
  new_order <- match(new_names_order, colnames(my_oncomat))
  return(new_order)
}







### Define colors for mutation types
mutation_colors <- c(Nonsense_Mutation="#ad7aff",Missense_Mutation="#377EB8",Frame_Shift_Del="#4DAF4A",
         In_Frame_Ins="#ff008c",Splice_Site="#FF7F00",Multi_Hit="#FFFF33",Frame_Shift_Ins="#A65628",
         In_Frame_Del="#f781bf",Translation_Start_Site="#400085",Nonstop_Mutation="#b68dfc",
         no_variants="#d6d6d6")
names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
### List defining functions for color and shape of cells in oncoplot
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # "0" = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
  #             gp = gpar(fill = "#CCCCCC", col = NA))
  # },
  "Nonsense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Nonsense Mutation"], col = NA))
  },
  "Missense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Missense Mutation"], col = NA))
  },
  "Frame Shift Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Frame Shift Del"], col = NA))
  },
  "In Frame Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["In Frame Ins"], col = NA))
  },
  "Splice Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Splice Site"], col = NA))
  },
  "Multi Hit" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Multi Hit"], col = NA))
  },
  "Frame Shift Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Frame Shift Ins"], col = NA))
  },
  "In Frame Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["In Frame Del"], col = NA))
  },
  "Nonstop Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Nonstop Mutation"], col = NA))
  },
  "Translation Start Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Translation Start Site"], col = NA))
  },
  "no variants" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              # gp = gpar(fill = "#e0e0e0", col = NA))
              gp = gpar(fill = "#CCCCCC", col = NA))
  }
)

