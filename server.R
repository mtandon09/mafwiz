library(shiny)
library(shinyjs)
library(openxlsx)
library(shinycssloaders)
library(shinyWidgets)
library(shinythemes)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(colourpicker)

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

options(shiny.maxRequestSize=300*1024^2) 

shinyServer(function(input, output, session) {
    
  plot_values <- reactiveValues()
  plot_values$raw_maf_obj <- NULL
  plot_values$maf_obj <- NULL
  plot_values$gene_interaction_pdf <- "maftools_somatic_interactions.png"
  # plot_values$oncoplot <- NULL
  # plot_values$oncoplot <- NULL
  
  observeEvent(input$get_tcga_data, {
    req(input$tcga_dataset)
    req(input$tcga_pipeline)
    
    save_folder="data/tcga_data"
    maf_file=paste0(save_folder,"/",input$tcga_dataset,".",input$tcga_pipeline,".maf")
    
    if (! file.exists(maf_file)) {
      if (!dir.exists(save_folder)) {dir.create(save_folder)}
      maf_download_progress <- shiny::Progress$new(session,max=100)
      maf_download_progress$set(value = 10, message = paste0("Downloading ",input$tcga_dataset,"..."))
      tcga_maf <- GDCquery_Maf(gsub("TCGA-","",input$tcga_dataset), 
                               pipelines = input$tcga_pipeline, 
                               directory = save_folder)
      write.table(tcga_maf, file=maf_file, quote=F, sep="\t", row.names = F, col.names = T)
      maf_download_progress$set(value = 90, message = paste0("Finished downloading ",input$tcga_dataset,"..."))
      maf_download_progress$close()
    }
    
    # plot_values$raw_maf_obj <- read_maf(session=session, maf_file)
    raw_maf <- read_maf(session=session, maf_file)
    
    plot_values$raw_maf_obj <- raw_maf
    
    click("run_apply_filters")
    # updateTabsetPanel(session, "main_tabs",
    #                   selected = "Visualizations")
  })
  
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
    } else {
      updateRadioButtons(session=session, inputId="burden_plot_type",selected = "Barplot")
    }
    plot_values$maf_obj <- filtered_maf
    
    updateTabsetPanel(session, "main_tabs",
                      selected = "Visualizations")
    
  })
  
  observeEvent(input$maf_file_upload, {
    
    req(input$maf_file_upload)
    myfile_path <- input$maf_file_upload$datapath
    # file_extension <- tools::file_ext(myfile_path)
    # plot_values$maf_obj <- read_maf(session=session, myfile_path)
    raw_maf <- read_maf(session=session, myfile_path)
    
    # raw_maf <- read_maf(session=session, maf_file)
    
    plot_values$raw_maf_obj <- raw_maf
    
    click("run_apply_filters")
    # updateTabsetPanel(session, "main_tabs", 
    #                   selected = "Visualizations")
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
      
      # ribbon_results <- make_single_ribbon_plot(plotopts$maf_obj, 
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
  
  output$oncoplot_output <- renderPlot({
    req(plot_values$maf_obj)
    oncoplot_progress <- shiny::Progress$new(session,max=100)
    oncoplot_progress$set(value = 5, message = "Making oncoplot...")
    myoncoplot <- make_oncoplot(plot_values$maf_obj, cohort_freq_thresh = input$onco_cohort_frac)
    # validate(
    #   need(class(oncoplot)=="ComplexHeatmap")
    # )
    oncoplot_progress$set(value = 100, message = "Returning plot...")
    oncoplot_progress$close()
    print(class(myoncoplot))
    if (class(myoncoplot)=="Heatmap") {
      draw(myoncoplot)
      plot_values$oncoplot <- myoncoplot
    }
  })  
  output$burden_output <- renderPlot({
    req(plot_values$maf_obj)
    burden_progress <- shiny::Progress$new(session,max=100)
    burden_progress$set(value = 5, message = "Making burden plot...")

    myburdenplot <- make_burden_plot(plot_values$maf_obj, input$burden_plot_type)
    burden_progress$set(value = 100, message = "Returning plot...")
    burden_progress$close()
    # print(class(myburdenplot[[2]]))
    # print(class(myburdenplot))
    plot_values$burdenplot <- myburdenplot
    # myburdenplot
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
                                               progress_func=updateProgress)
    mutsig_progress$set(value = 100, message = "Returning plot...")
    mutsig_progress$close()
    print(class(mymutsigplot))
    plot_values$mutsigplot <- mymutsigplot
    draw(mymutsigplot)
  })
  
  
  output$generibbon_output <- renderPlot({
    req(plot_values$maf_obj)
    # req(input$gene_comut_topN)s
    # print(input$gene_comut_customribboncolor)
    # print(input$gene_comut_ribbon_color)
    # mycoloropt <- NULL
    # if (input$gene_comut_customribboncolor) {
    #   mycoloropt <- input$gene_comut_ribbon_color
    # }
    
    # print(isolate(input$gene_comut_topN))
    
    # plotopts <- list(maf_obj=plot_values$maf_obj, 
    #                  ribbon_color=mycoloropt, 
    #                  topN=plot_values$gene_comut_topN,
    #                  pval_low=plot_values$gene_comut_pvalLow, 
    #                  pval_high=plot_values$gene_comut_pvalHigh,
    #                  plot_file=plot_values$gene_interaction_pdf,
    #                  plot_frac_mut_axis=TRUE,
    #                  scale_ribbon_to_fracmut=TRUE)
    
    
    # mycoloropt <- NULL
    # if (input$gene_comut_customribboncolor) {
    #   mycoloropt <- input$gene_comut_ribbon_color
    # }
    # 
    # plotopts <- list(maf_obj=plot_values$maf_obj, 
    #                  ribbon_color=mycoloropt, 
    #                  topN=isolate(input$gene_comut_topN),
    #                  pval_low=isolate(input$gene_comut_pvalLow), 
    #                  pval_high=isolate(input$gene_comut_pvalHigh),
    #                  plot_file=isolate(plot_values$gene_interaction_pdf),
    #                  plot_frac_mut_axis=TRUE,
    #                  scale_ribbon_to_fracmut=TRUE)
    # 
    # ribbon_progress <- shiny::Progress$new(session,max=100)
    # ribbon_progress$set(value = 5, message = "Making co-mutated gene ribbon plot...")
    # 
    # updateProgress <- function(value = NULL, detail = NULL) {
    #   if (is.null(value)) {
    #     value <- ribbon_progress$getValue()
    #     value <- value + (ribbon_progress$getMax() - value) / 5
    #   }
    #   ribbon_progress$set(value = value, message = detail)
    # }
    # 
    # # ribbon_results <- make_single_ribbon_plot(plotopts$maf_obj, 
    # make_single_ribbon_plot(plotopts$maf_obj, 
    #                                           ribbon_color=plotopts$ribbon_color, 
    #                                           topN=plotopts$topN,
    #                                           pval_low=plotopts$pval_low, 
    #                                           pval_high=plotopts$pval_high,
    #                                           plot_file=plotopts$plot_file,
    #                                           progress_func=plotopts$progress_func,
    #                                           plot_frac_mut_axis=plotopts$plot_frac_mut_axis,
    #                                           scale_ribbon_to_fracmut=plotopts$scale_ribbon_to_fracmut)  
    # 
    # ribbon_progress$set(value = 100, message = "Rendering plot...")
    # ribbon_progress$close()
    
    # draw(ribbon_results[[1]])
    # draw(ribbon_results[[2]], x = unit(0.5, "npc"), y = unit(0.05, "npc"), just = c("center"))
    
    
    click("plot_generibbon")
    # make_gene_ribbon_plot(session, plotopts)
    # print(class(mymutsigplot))
    # plot_values$mutsigplot <- mymutsigplot
    # draw(mymutsigplot)
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
  
  observeEvent(input$clear_maf_file, {
    reset("maf_file_upload")
    plot_values$maf_obj <- NULL
  })
  
  observeEvent(input$load_example_data, {
    reset("maf_file_upload")
    myfile_path <- "data/TCGA-CHOL.maf"
    # plot_values$maf_obj <- read_maf(session=session, myfile_path)
    raw_maf <- read_maf(session=session, myfile_path)
    
    plot_values$raw_maf_obj <- raw_maf
    
    click("run_apply_filters")
    # updateTabsetPanel(session, "main_tabs",
    #                   selected = "Visualizations")
  })
  
    
    
  }
)