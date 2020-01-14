library(shiny)
library(shinyjs)
library(shinythemes)
library(ggbeeswarm)
library(colourpicker)

library(shinyBS)
library(shinycssloaders)
library(shinyWidgets)
library(shinythemes)

source("helper_functions.oncoplot.R")
source("helper_functions.shiny.R")

# library(DT)
# dataset <- diamonds
# Define UI for data upload app ----
shinyUI(
  navbarPage(theme = shinytheme("flatly"),windowTitle="MafViz", selected="Upload and Filter",
             title="MafWiz", id="main_tabs",
             useShinyjs(),
             tabPanel("Upload and Filter",
                      h1("Select a MAF file"),
                      tags$hr(),
                      fluidRow(
                        column(6, #offset = 1,
                               h4("Use TCGA data"),
                               selectizeInput("tcga_dataset", "TCGA Dataset",
                                           choices=TCGA_project_choices),
                               # tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 20px;}"),
                               selectInput("tcga_pipeline", "Variant Calling Pipeline",
                                            choices=TCGA_pipeline_choices, selected="mutect2",
                                            selectize = F),
                               actionButton("get_tcga_data","Download"),
                        ),
                        column(4,
                               h4("Use your own data"),
                              actionLink("clear_maf_file", "Reset",
                                         style="color: #cccccc; font-size: 10spx;"),
                              fileInput("maf_file_upload", "Upload MAF file",
                                        multiple = FALSE,
                                        accept = c("text",".maf",".txt")),
                              actionLink("load_example_data", "Use example data (TCGA-CHOL)",
                                         style="color: #0000ff; font-size: 12spx;"),
                               
                        )
                      ),
                      tags$br(),
                      tags$br(),
                      tags$hr(),
                      h3("Apply filters"),
                      # disabled(
                        # checkboxInput("apply_filters", "Filter MAF file?", FALSE),
                        checkboxInput("exclude_flags", "Exclude known frequently mutated genes (FLAGs)?", FALSE),
                        actionLink("flags_info","More info on FLAGs genes"),
                        actionButton("run_apply_filters","Apply Filters"),
                      # ),
                      tags$hr(),
                      h3("Add Sample Information"),
                      disabled(
                        fileInput("sample_file_upload", "Choose file",
                                  multiple = FALSE,
                                  accept = c("text",".txt",".xlsx",".xls",".csv"))
                      )
              ),
             
             tabPanel("Visualizations",
                      tabsetPanel(id = "viz_panel", type="pills",
                                  tabPanel("Mutation Burden",
                                           tags$hr(),
                                           sidebarPanel("", width=3,
                                              radioButtons("burden_plot_type",label = "Type of plot", 
                                                           choices = c("Barplot", "Dotplot"))
                                           ),
                                           mainPanel(
                                             withSpinner(plotOutput("burden_output", width = "90%", height = "600px", click = NULL,
                                                                    dblclick = NULL, hover = NULL, hoverDelay = NULL,
                                                                    hoverDelayType = NULL, brush = NULL, clickId = NULL,
                                                                    hoverId = NULL, inline = FALSE), type=1),
                                             tags$hr(),
                                             downloadButton("download_burdenplot_button",label="Download plot"),
                                             bsModal("download_burdenplot_modal","Download plot","download_burdenplot_button",
                                                     radioButtons("burdenplot_save_type","Format",c("pdf","png","tiff","svg"),selected="pdf"),
                                                     numericInput("burdenplot_save_width", "Width (in)",value=8, min=2, max=50, step=0.5),
                                                     numericInput("burdenplot_save_height", "Height (in)",value=8, min=2, max=50, step=0.5),
                                                     downloadButton("download_burdenplot","Download"))
                                           )
                                  ),
                                  tabPanel("Oncoplot",
                                           tags$hr(),
                                           # sidebarPanel("Plot Options", width=2),
                                           mainPanel(
                                             tags$hr(),
                                             # withSpinner(plotOutput("oncoplot_output"),type=1),
                                             tabsetPanel(id = "viz_panel", type="tabs",
                                              tabPanel("Plot",
                                                       withSpinner(plotOutput("oncoplot_output"),type=1)
                                              ),
                                              tabPanel("Options",
                                                       # p("Plot Options"),
                                                       numericInput("onco_cohort_frac","Plot genes with mutation in fraction of cohort",
                                                                    value=0.05, min=0, max=1, step=0.01)
                                              )
                                             ),
                                              
                                             tags$hr(),
                                             downloadButton("download_oncoplot_button",label="Download plot"),
                                             bsModal("download_oncoplot_modal","Download plot","download_oncoplot_button",
                                                     radioButtons("oncoplot_save_type","Format",c("pdf"),selected="pdf"),
                                                     numericInput("oncoplot_save_width", "Width (in)",value=8, min=2, max=50, step=0.5),
                                                     numericInput("oncoplot_save_height", "Height (in)",value=8, min=2, max=50, step=0.5),
                                                     downloadButton("download_oncoplot","Download"))
                                           )
                                  ),
                                  tabPanel("Mutational Signatures",
                                           sidebarPanel("", width=3,
                                                        # selectInput("genome_select",label="Genome Build",choices=c(NULL,"hg19","hg38")),selectize=F,
                                                        checkboxInput("use_syn_mut", "Use Synonymous Mutations?", FALSE)
                                           ),
                                           mainPanel(
                                             withSpinner(plotOutput("mutsig_output", width = "90%", height = "600px", click = NULL,
                                                                    dblclick = NULL, hover = NULL, hoverDelay = NULL,
                                                                    hoverDelayType = NULL, brush = NULL, clickId = NULL,
                                                                    hoverId = NULL, inline = FALSE), type=1),
                                             tags$hr(),
                                             downloadButton("download_mutsigplot_button",label="Download plot"),
                                             bsModal("download_mutsigplot_modal","Download plot","download_mutsigplot_button",
                                                     radioButtons("mutsigplot_save_type","Format",c("pdf"),selected="pdf"),
                                                     numericInput("mutsigplot_save_width", "Width (in)",value=8, min=2, max=50, step=0.5),
                                                     numericInput("mutsigplot_save_height", "Height (in)",value=8, min=2, max=50, step=0.5),
                                                     downloadButton("download_mutsigplot","Download"))
                                           )
                                  ),
                                  tabPanel("Co-mutated Genes",
                                           sidebarPanel("", width=3,
                                                        # selectInput("genome_select",label="Genome Build",choices=c(NULL,"hg19","hg38")),selectize=F,
                                                        checkboxInput("use_syn_mut_genecomut", "Use Synonymous Mutations?", FALSE),
                                                        bsCollapse(id = "gene_comut_options", open = NULL,
                                                                   bsCollapsePanel(title="Choose Genes", 
                                                                                   style="info",
                                                                                   actionLink("reset_gene_comut_options", "Reset",
                                                                                              # style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                                                                              # style="color: #333333; font-weight: bold; font-size: 20px; background-color: #ffcc00; border-color: #2e6da4"),
                                                                                              # style="color: #000000; font-size: 6px; background-color: #009900; border-color: #2e6da4"),
                                                                                              style="color: #cccccc; font-size: 10spx;"),
                                                                                   # Input: Select a file ----
                                                                                   
                                                                                   numericInput("gene_comut_topN","Analyze this many top mutated genes", value=25, min=3, max=100, step=1)
                                                                                   
                                                                   ),
                                                                   bsCollapsePanel(title="Plot Options",
                                                                                   style="info",
                                                                                   numericInput("gene_comut_pvalHigh","Minimum p-value of co-occurence for plot", value=0.1, min=0, max=1, step=1),
                                                                                   numericInput("gene_comut_pvalLow","p-value of co-occurence shown with dotted lines", value=0.05, min=0, max=1, step=1),
                                                                                   checkboxInput("gene_comut_customribboncolor", "Use custom ribbon color",value = FALSE),
                                                                                   colourInput("gene_comut_ribbon_color", "Select color for ribbon",value="steelblue",returnName=TRUE)
                                                                   )
                                                       ),
                                                       actionButton("plot_generibbon","Make plot")
                                           ),
                                           mainPanel(
                                             tabsetPanel(id = "generibbon_panel", type="tabs",
                                                         tabPanel("Ribbon Plot",
                                                                  withSpinner(plotOutput("generibbon_output", width = "90%", height = "600px", click = NULL,
                                                                                         dblclick = NULL, hover = NULL, hoverDelay = NULL,
                                                                                         hoverDelayType = NULL, brush = NULL, clickId = NULL,
                                                                                         hoverId = NULL, inline = FALSE), type=1)
                                                         ),
                                                         tabPanel("Interaction Matrix",
                                                                  # p("Plot Options"),
                                                                  withSpinner(imageOutput("genematrix_output", width = "90%", height = "600px", click = NULL,
                                                                                         dblclick = NULL, hover = NULL, hoverDelay = NULL,
                                                                                         hoverDelayType = NULL, brush = NULL, clickId = NULL,
                                                                                         hoverId = NULL, inline = FALSE), type=1)
                                                         )
                                             ),
                                             
                                             tags$hr(),
                                             downloadButton("download_generibbon_button",label="Download plot"),
                                             bsModal("download_generibbon_modal","Download plot","download_generibbon_button",
                                                     radioButtons("generibbon_save_type","Format",c("pdf"),selected="pdf"),
                                                     numericInput("generibbon_save_width", "Width (in)",value=8, min=2, max=50, step=0.5),
                                                     numericInput("generibbon_save_height", "Height (in)",value=8, min=2, max=50, step=0.5),
                                                     downloadButton("download_generibbon_plot","Download"))
                                           )
                                  )
                      )
             ),
             tabPanel("Table Output",
                      p("Nothing here yet")
             ),
             navbarMenu("More",
                        tabPanel("Stuff to add in the future",
                                 p("-- Launch MutSigCV?")))
  )
)