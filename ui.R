library(shiny)
library(shinyjs)
library(shinyBS)
library(shinycssloaders)
library(shinyWidgets)
library(shinythemes)
library(DT)
library(rhandsontable)
library(colourpicker)
# library(shinybusy)
library(shinydashboard)

source("helper_functions.oncoplot.R")
source("helper_functions.shiny.R")

shinyUI(
  
  dashboardPage(skin = "black",
                
      dashboardHeader(title = "MafWiz",titleWidth=100),
      dashboardSidebar(
        sidebarMenu(
          # Setting id makes input$tabs give the tabName of currently-selected tab
          id = "tab_menu",
          menuItem("Get Data", icon=icon("cloud-download"), tabName = "get-data", selected=T),
          menuItem("Pick Clinical Variables", tabName = "get-clinical-variables", icon = icon("search")),
          menuItem("Visualizations", icon = icon("bar-chart-o"), tabName = "visualizations"),
          menuItem("Table Output", icon = icon("table"), tabName = "table-output"),
          menuItem("VCF2MAF", icon = icon("random"), tabName = "vcf2maf", badgeLabel = "new",
                   badgeColor = "green")
        )
      ),
      dashboardBody(
             tabItems(
                 # tabPanel("Upload and Filter",
                 tabItem(tabName="get-data",
                         useShinyjs(),
                          h1("Select a MAF file"),
                          # tags$hr(),
                          fluidRow(
                            column(6, #offset = 1,
                                   h4("Use TCGA data"),
                                   selectizeInput("tcga_dataset", "TCGA Dataset",
                                               choices=TCGA_project_choices),
                                   # tags$style(type='text/css', ".selectize-input { font-size: 12px; line-height: 20px;}"),
                                   selectInput("tcga_pipeline", "Variant Calling Pipeline",
                                                choices=TCGA_pipeline_choices, selected="mutect2",
                                                selectize = F),
                                   actionButton("get_tcga_data","Download", icon=icon("cloud-download")),
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
                          h3("Add Sample Information"),
                          # h3(actionBttn("pick_clin_vars", label="Add Sample Information", style="unite")),
                          # disabled(
                          fluidRow(
                            column(6, #offset = 1,
                                   h4("Get TCGA clinical data"),
                                   actionButton("get_tcga_clinical_data","Download", icon=icon("cloud-download")),
                            ),
                            column(4,
                                   h4("Use your own data"),
                                   actionLink("clear_clinical_data_file", "Reset",
                                              style="color: #cccccc; font-size: 10spx;"),
                                   fileInput("sample_file_upload", "Choose Tab-separated text or Excel file",
                                            multiple = FALSE,
                                            accept = c("text",".txt",".xlsx",".xls",".csv")),
                                   actionLink("load_example_clinical_data", "Use example data (TCGA-CHOL)",
                                              style="color: #0000ff; font-size: 12spx;"),
                                   
                            )
                          ),
                          tags$br(),
                          tags$hr(),
                          h3("Apply filters"),
                          checkboxInput("exclude_flags", "Exclude known frequently mutated genes (FLAGs)?", TRUE),
                          actionLink("flags_info","More info on FLAGs genes"),
                          bsModal("flags_info_popup","FLAGs filter","flags_info",
                                  # hr(),
                                  includeMarkdown("html_content/flags.md"),
                                  size="large"),
                          actionButton("run_apply_filters","Apply Filters"),
                          tags$hr(),
                          # actionButton("pick_clin_vars","Pick Clinical Variables"),
                          tags$hr()
                          
                          # )
                  ),
                 tabItem(tabName="get-clinical-variables",
                          h1("Sample Information Table"),
                          h2(actionBttn("back","Go to Visualizations", icon=icon("chart-bar"))),
                          tags$hr(),
                          fluidRow(
                            column(3, #offset = 1,
                                   h4("Find clinical variable"),
                                   p("try searching for age, gender, stage, etc."),
                                   # actionButton("add_curr_var_data","Add to sample table"),
                                   actionBttn("add_curr_var_data","Add to sample table", icon=icon("plus-square"), style="fill", size="sm"),
                                   tags$hr(),
                                   selectInput("select_curr_clin_var", "Column Name",
                                               choices=c(), selected=NULL,
                                               selectize = T),
                                   tags$hr(),
                                   DTOutput("curr_clin_var_table", height = "300px")
                            ),
                            column(8,
                                   h4("Variables for analysis"),
                                   actionLink("clear_clinical_data", "Reset",
                                              style="color: #cccccc; font-size: 10spx;"),
                                   tags$style(HTML('table.dataTable tr:nth-child(even) {background-color: #deebff !important;}')),
                                   tags$style(HTML('table.dataTable tr:nth-child(odd) {background-color: #fafcff !important;}')),
                                   tags$style(HTML('table.dataTable th {background-color: white !important;}')),
                                   # dataTableOutput("clin_dat_header_table"),
                                   # dataTableOutput("maf_clin_dat_table")
                                   rHandsontableOutput("maf_clin_dat_table", height = "400px")
                            )
                          ),
                          tags$hr(),
                          h2("Annotation Legend"),
                          actionBttn("pick_anno_colors","Pick Colors", icon=icon("palette"), style="jelly", color="default"),
                          bsModal("anno_color_picker","Pick Annotation Colors","pick_anno_colors",size="large",
                                  selectInput("color_var_picker","Select Variable", choices = c()),
                                  selectInput("color_var_type_picker","Variable Type", choices = c("Category","Numeric")),
                                  textOutput("test_txtout"),
                                  conditionalPanel(
                                    condition = "input.color_var_type_picker == 'Numeric'",
                                    fixedRow(column(1,colourInput("curr_min_color", "Min", "white",showColour="text")),
                                             column(1,colourInput("curr_max_color", "Max", "Black",showColour="text"))),
                                  ),
                                  plotOutput("curr_anno_legend",width = "90%", height = "300px")
                                  ),
                          withSpinner(plotOutput("anno_legend_preview", width = "80%", height = "100px"),type=1)
                          
                 ),
                 tabItem(tabName="vcf2maf",
                         h2("Convert VCF to MAF")
                 ),
                 tabItem(tabName="visualizations",
                          tabsetPanel(id = "viz_panel", type="tabs",
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
                                                           withSpinner(plotOutput("oncoplot_output", width = "800px", height = "600px"),type=1)
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
                                                 withSpinner(plotOutput("mutsig_output", width = "800px", height = "800px", click = NULL,
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
                 tabItem(tabName="table-output",
                          # p("Nothing here yet")
                          h1("Variant Table"),
                          tags$hr(),
                          # withSpinner(rHandsontableOutput("output_variant_table_rhandson"),type=1)
                          withSpinner(DTOutput("output_variant_table_DT"),type=1)
                 )
        )
      )
  )
)