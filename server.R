library(shiny)
library(shinyjs)
library(openxlsx)
library(shinycssloaders)
library(shinyWidgets)
library(shinythemes)
library(RColorBrewer)
library(ggplot2)
library(reshape2)

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
  plot_values$maf_obj <- NULL
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
    
    plot_values$maf_obj <- read_maf(session=session, maf_file)
    updateTabsetPanel(session, "main_tabs",
                      selected = "Visualizations")
  })
  
  observeEvent(input$maf_file_upload, {
    
    req(input$maf_file_upload)
    myfile_path <- input$maf_file_upload$datapath
    # file_extension <- tools::file_ext(myfile_path)
    plot_values$maf_obj <- read_maf(session=session, myfile_path)
    
    updateTabsetPanel(session, "main_tabs", 
                      selected = "Visualizations")
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
    plot_values$maf_obj <- read_maf(session=session, myfile_path)
    updateTabsetPanel(session, "main_tabs",
                      selected = "Visualizations")
  })
  
    
    
  }
)