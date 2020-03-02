library(shiny)
library(shinyjs)
library(DT)
library(openxlsx)
library(shinycssloaders)
library(shinyWidgets)
library(shinythemes)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(colourpicker)
library(rhandsontable)
library(shinybusy)

library(maftools)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(ggbeeswarm)
library(TCGAbiolinks)
library(NMF)
library(MutationalPatterns)

library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)

source("helper_functions.oncoplot.R")
source("helper_functions.shiny.R")

options(shiny.maxRequestSize=1000*1024^2)  ## Increase upload limit to 1G 

shinyServer(function(input, output, session) {
    
  plot_values <- reactiveValues()
  plot_values$raw_maf_obj <- NULL
  plot_values$maf_obj <- NULL
  plot_values$gene_interaction_pdf <- "maftools_somatic_interactions.png"
  plot_values$output_table <- NULL
  
  clin_var_values <- reactiveValues()
  clin_var_values$clin_var_data <- NULL
  clin_var_values$clin_data_upload <- NULL
  clin_var_values$clin_data_file <- NULL
  clin_var_values$clin_data_columns <- NULL
  clin_var_values$clin_anno_colors <- NULL
  
  observeEvent(input$get_tcga_data, {
    req(input$tcga_dataset)
    req(input$tcga_pipeline)
    
    save_folder="data/tcga_data"
    maf_file=file.path(save_folder,input$tcga_dataset,paste0(input$tcga_dataset,".",input$tcga_pipeline,".maf"))
    
    if (! file.exists(maf_file)) {
      if (!dir.exists(dirname(maf_file))) {dir.create(dirname(maf_file), recursive = T)}
      maf_download_progress <- shiny::Progress$new(session,max=100)
      maf_download_progress$set(value = 10, message = paste0("Downloading ",input$tcga_dataset,"..."))
      tcga_maf <- GDCquery_Maf(gsub("TCGA-","",input$tcga_dataset), 
                               pipelines = input$tcga_pipeline, 
                               directory = save_folder)
      tcga_maf$Tumor_Sample_Barcode_original <- tcga_maf$Tumor_Sample_Barcode
      tcga_maf$Tumor_Sample_Barcode <-unlist(lapply(strsplit(tcga_maf$Tumor_Sample_Barcode, "-"), function(x) {paste0(x[1:3], collapse="-")}))
      write.table(tcga_maf, file=maf_file, quote=F, sep="\t", row.names = F, col.names = T)
      maf_download_progress$set(value = 90, message = paste0("Finished downloading ",input$tcga_dataset,"..."))
      maf_download_progress$close()
    }
    
    raw_maf <- read_maf(session=session, maf_file)
    
    plot_values$raw_maf_obj <- raw_maf
    click("get_tcga_clinical_data")
    click("run_apply_filters")
  })
  
  observeEvent(input$get_tcga_clinical_data, {
    req(input$tcga_dataset)
    print("get tcga clin data")
    click("clear_clinical_data")
    
    save_folder="data/tcga_data"
    tcga_clinical_file=file.path(save_folder,input$tcga_dataset,paste0(input$tcga_dataset,".clinical.txt"))
    if (! file.exists(tcga_clinical_file)) {
      if (!dir.exists(dirname(tcga_clinical_file))) {dir.create(dirname(tcga_clinical_file), recursive = T)}
      tcga_clinical <- GDCquery_clinic(project = input$tcga_dataset, type = "clinical")
      write.table(tcga_clinical, file=tcga_clinical_file, quote=T, sep="\t", row.names = F, col.names = T)
    }
    clin_var_values$clin_data_file <- tcga_clinical_file
    
  })
    
  observeEvent(input$sample_file_upload, {
    req(input$tcga_dataset)
    print("get local clin data")
    click("clear_clinical_data")
    clin_var_values$clin_data_file <- input$sample_file_upload$datapath
    
  })
  
  
  observeEvent(input$load_example_clinical_data, {
    reset("sample_file_upload")
    click("clear_clinical_data")
    myfile_path <- file.path("data","TCGA-CHOL.clinical.txt")
    
    print("get example clin data")
    clin_var_values$clin_data_file <- myfile_path
    
  })
  
  
  
  observe({
    req(clin_var_values$clin_data_file)
    click("clear_clinical_data")
    clinical_file <- clin_var_values$clin_data_file 
    file_type=tools::file_ext(clinical_file)
    if (grepl("xls*",file_type)) {
      clin_data_upload <- read.xlsx(clinical_file)
    } else if (grepl("txt|tsv", file_type)) {
      clin_data_upload <- read.table(clinical_file, sep="\t", header=T)
    } else {
      stop(paste0("Don't know what to do with file extension: ", file_type))
    }
    
    clin_var_values$clin_data_upload <- clin_data_upload
    
    updateSelectInput(session, "select_curr_clin_var", choices = colnames(clin_var_values$clin_data_upload))
    
    updateTabItems(session, "tab_menu",
                      selected = "get-clinical-variables")

    
  })
  
  
  observeEvent(input$add_curr_var_data, {
    print("adding new var to clin data")
    
    curr_data <- clin_var_values$clin_data_upload
    # curr_data <- input$maf_clin_dat_table_cell_info
    print(input$maf_clin_dat_table_cell_info)
    curr_col <- match(input$select_curr_clin_var,colnames(clin_var_values$clin_data_upload))
    curr_data_table <- curr_data[,curr_col,drop=F]
    
    if(is.null(clin_var_values$clin_var_data)){
      clin_var_values$clin_var_data <- curr_data_table
    } else {
      clin_var_values$clin_var_data <- cbind(clin_var_values$clin_var_data,curr_data_table)
    }
    updateSelectInput(session,"color_var_picker","Select Variable",choices=colnames(clin_var_values$clin_var_data))
    
    
  })
  
  observeEvent(input$clear_clinical_data, {
    clin_var_values$clin_var_data <- NULL
    clin_var_values$clin_anno_colors <- NULL
  })
  
  
  observeEvent(input$color_var_picker, {
    req(clin_var_values$clin_anno_colors)
    print(clin_var_values$clin_anno_colors[[input$color_var_picker]])
    curr_clin_var_color <- clin_var_values$clin_anno_colors[[input$color_var_picker]]
    curr_type <- ifelse(is.function(curr_clin_var_color),"Numeric","Category")
    updateSelectInput(session,inputId = "color_var_type_picker",selected=curr_type)
    
  })
  
  output$test_txtout <- renderText(input$color_var_picker)

  
  
  
  
  
  
  
  observeEvent(input$run_apply_filters, {
    req(plot_values$raw_maf_obj)
    raw_maf <- plot_values$raw_maf_obj
    
    if (input$exclude_flags) {
      flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
      genes_to_keep <- setdiff(raw_maf@gene.summary$Hugo_Symbol, flag_genes)
      filtered_maf <- subsetMaf(raw_maf,genes=genes_to_keep, mafObj = T)
    } else {
      filtered_maf <- raw_maf
    }
    
    if (nrow(filtered_maf@gene.summary) > 20) {
      print("here")
      updateRadioButtons(session=session, inputId="burden_plot_type",selected = "Dotplot")
    }
    plot_values$maf_obj <- filtered_maf
    
    plot_values$output_table <- make_variant_table(plot_values$maf_obj)
    
    updateTabItems(session, "tab_menu",
                      selected = "visualizations")

    
  })
  
  observeEvent(input$maf_file_upload, {
    
    req(input$maf_file_upload)
    myfile_path <- input$maf_file_upload$datapath
    raw_maf <- read_maf(session=session, myfile_path)
    
    plot_values$raw_maf_obj <- raw_maf
    
    click("run_apply_filters")
  })
  
  observeEvent(input$plot_generibbon, {
    output$generibbon_output <- renderPlot({
      req(plot_values$maf_obj)
      
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
      
      ribbon_progress$set(value = 100, message = "Rendering plot...")
      ribbon_progress$close()
    })
    
  })
  
  observeEvent(input$gene_comut_customribboncolor, {
    
    if (isolate(input$gene_comut_customribboncolor)) {
      enable(input$gene_comut_ribbon_color)
    } else {
      disable(input$gene_comut_ribbon_color)
    }
    
  })
  
  
  
  
  
  
  
  
  
  
  
  
  output$curr_clin_var_table <- renderDT({
    req(clin_var_values$clin_data_upload)
    curr_data <- clin_var_values$clin_data_upload
    curr_col <- match(input$select_curr_clin_var,colnames(clin_var_values$clin_data_upload))
    curr_data_table <- curr_data[,curr_col,drop=F]
    datatable(curr_data_table,
              rownames=F,
              editable=F,
              extensions = c("AutoFill","ColReorder","KeyTable","Scroller"),
              options = list(searching = T,
                             deferRender = T,
                             autoFill = F,
                             colReorder = F,
                             keys = F,
                             scrollY = 200,
                             scroller = T)
              ) %>% formatStyle(names(curr_data_table), backgroundColor = "#fafcff")
  },
  server=F
  )
  
  
  output$maf_clin_dat_table <- renderRHandsontable({
    req(clin_var_values$clin_var_data)
    curr_data <- clin_var_values$clin_var_data
    curr_data <- curr_data[,!duplicated(colnames(curr_data)), drop=F]
    
    rhandsontable(curr_data) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE)
  })

  
  output$output_variant_table_rhandson <- renderRHandsontable({
    req(plot_values$output_table)
    
    
    var_table <- plot_values$output_table
    
    rhandsontable(var_table)
  })
  
  output$output_variant_table_DT <- renderDT({
    req(plot_values$output_table)
    
    
    var_table <- plot_values$output_table
    
    datatable(var_table,
              extensions = c("AutoFill","ColReorder","KeyTable","Scroller"),
              options = list(searching = T,
                             deferRender = T,
                             autoFill = F,
                             colReorder = F,
                             keys = F,
                             scrollY = 400,
                             scroller = T)
    )
  },
    server=F
  )
  
  output$anno_legend_preview <- renderPlot({
    req(clin_var_values$clin_var_data)
    
    # print(colnames(clin_var_values$clin_var_data))
    
    anno_data <- clin_var_values$clin_var_data
    make_anno <- apply(anno_data,2,function(x) {
      ret_val=TRUE
      if (is.factor(x) && levels(x)>10) {
        ret_val=FALSE 
      }
      return(ret_val)
    })
    anno_data <- anno_data[,make_anno, drop=F]
    if (ncol(anno_data) > 0) {
      if(is.null(clin_var_values$clin_anno_colors)) {
        myanno <- HeatmapAnnotation(df=anno_data,
                                    which="column",
                                    show_annotation_name = TRUE,
                                    annotation_name_side = "right")
        
      } else {
        myanno <- HeatmapAnnotation(df=clin_var_values$clin_var_data,
                                  which="column",
                                  col = clin_var_values$clin_anno_colors,
                                  show_annotation_name = TRUE,
                                  annotation_name_side = "right")
      }
      
      mycolors <- lapply(names(myanno@anno_list), function(x) {
                          if (myanno@anno_list[[x]]@color_mapping@type=="continuous") {
                            ret_val <- myanno@anno_list[[x]]@color_mapping@col_fun
                          } else {
                            ret_val <- myanno@anno_list[[x]]@color_mapping@colors
                          }
                          
                                       
                          return(ret_val)
      })
      names(mycolors) <- names(myanno@anno_list)
      clin_var_values$clin_anno_colors <- mycolors
      draw(myanno)
    }
    
    
  }) 
  
  
  output$oncoplot_output <- renderPlot({
    req(plot_values$maf_obj)
    oncoplot_progress <- shiny::Progress$new(session,max=100)
    oncoplot_progress$set(value = 5, message = "Making oncoplot...")
    myoncoplot <- make_oncoplot(plot_values$maf_obj, cohort_freq_thresh = input$onco_cohort_frac,
                                clin_data = clin_var_values$clin_var_data,
                                clin_data_colors = clin_var_values$clin_anno_colors
                                )

    oncoplot_progress$set(value = 100, message = "Returning plot...")
    if (grepl("Heatmap",class(myoncoplot))) {
      draw(myoncoplot)
      plot_values$oncoplot <- myoncoplot
    }
    oncoplot_progress$close()
  })  
  
  output$burden_output <- renderPlot({
    req(plot_values$maf_obj)
    burden_progress <- shiny::Progress$new(session,max=100)
    burden_progress$set(value = 5, message = "Making burden plot...")

    myburdenplot <- make_burden_plot(plot_values$maf_obj, input$burden_plot_type)
    burden_progress$set(value = 100, message = "Returning plot...")
    burden_progress$close()

    plot_values$burdenplot <- myburdenplot
    plot(myburdenplot)
  })
  
  output$mutsig_output <- renderPlot({
    req(plot_values$maf_obj)
    mutsig_progress <- shiny::Progress$new(session,max=100)
    mutsig_progress$set(value = 5, message = "Making mutational signatures plot...")
    
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- mutsig_progress$getValue()
        value <- value + (mutsig_progress$getMax() - value) / 5
      }
      mutsig_progress$set(value = value, message = detail)
    }
    # print(input$genome_select)
    mymutsigplot <- make_mut_signature_heatmap(plot_values$maf_obj, 
                                               # genome_build=input$genome_select,
                                               clin_data = clin_var_values$clin_var_data,
                                               clin_data_colors = clin_var_values$clin_anno_colors,
                                               progress_func=updateProgress)
    mutsig_progress$set(value = 100, message = "Returning plot...")
    mutsig_progress$close()
    print(class(mymutsigplot))
    plot_values$mutsigplot <- mymutsigplot
    draw(mymutsigplot)
  })
  
  output$generibbon_output <- renderPlot({
    req(plot_values$maf_obj)
    
    click("plot_generibbon")
  })
  
  output$genematrix_output <- renderImage({
    validate(need(file.exists(plot_values$gene_interaction_pdf),"No plot file found."))
    print(plot_values$gene_interaction_pdf)
    width  <- session$clientData$output_image_width
    height <- session$clientData$output_image_height
    list(
      src = plot_values$gene_interaction_pdf,
      contentType = "image/png",
      width = width,
      height = height,
      alt = "Output from maftools' somaticInteractions()"
    )
  },deleteFile=F)
  
  output$download_oncoplot <- downloadHandler (
    filename = function(){
      paste("oncoplot",input$oncoplot_save_type, sep=".")
    },
    content = function(ff) {
      req(plot_values$oncoplot)
      pdf(ff, height=input$oncoplot_save_height,width=input$oncoplot_save_width)
      draw(plot_values$oncoplot)
      dev.off()
    }
  )
  
  
  observeEvent(input$clear_maf_file, {
    reset("maf_file_upload")
    plot_values$maf_obj <- NULL
  })
  
  observeEvent(input$load_example_data, {
    reset("maf_file_upload")
    myfile_path <- file.path("data","TCGA-CHOL.maf")
    # plot_values$maf_obj <- read_maf(session=session, myfile_path)
    raw_maf <- read_maf(session=session, myfile_path)
    
    plot_values$raw_maf_obj <- raw_maf
    
    click("run_apply_filters")
  })

  
  
  output$download_burdenplot <- downloadHandler (
    filename = function(){
      paste("burden",input$burdenplot_save_type, sep=".")
    },
    content = function(ff) {
      req(plot_values$burdenplot)
      ggsave(plot_values$burdenplot, filename = ff,height=input$burdenplot_save_height,width=input$burdenplot_save_width,units="in")
      
    }
  )
  
  output$download_mutsigplot <- downloadHandler (
    filename = function(){
      paste("burden",input$mutsigplot_save_type, sep=".")
    },
    content = function(ff) {
      req(plot_values$mutsigplot)
      pdf(ff, height=input$mutsigplot_save_height,width=input$mutsigplot_save_width)
      draw(plot_values$mutsigplot)
      dev.off()
    }
  ) 
    
  }
)