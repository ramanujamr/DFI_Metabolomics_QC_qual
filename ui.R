# QUAL UI

ui <- fluidPage(shinytheme("journal"),
                useShinyjs(),
                titlePanel("DFI Metabolomics QC - Qual"),
                
                fluidRow(
                  h3("1. UPLOAD DATA"),
                  hr(),
                  column(width=3,
                         selectInput("filename", label="Select csv file to upload",
                                     list.files(wddir, pattern = "*.csv", ignore.case = T)),
                         actionButton("Button_refresh_csv", "Refresh", icon("sync"), width="100px"),
                         actionButton("Button_upload_csv", "Upload", icon("upload"), width="100px",
                                      style="color: #fff; background-color: #2346b0; border-color: #2e6da4"),
                         textOutput("Textout_filename"),
                         textOutput("Textout_panel")),
                  
                  fluidRow(
                    br(), br(),
                    column(width=3, align="center", offset=9,
                         shiny::actionButton("Button_generate_boxplots", "Generate Boxplots", icon("chart-bar"), width="200px",
                                             style="color: #fff; background-color: #2346b0; border-color: #2e6da4"))
                    )
                ),
                  br(), hr(),
                  
                  fluidRow(
                    h3("2. RAW DATA"),
                    br(),
                    column(width=12, align="center", plotOutput("Plot_boxplots", height="auto"))
                  ),
                
                br(), br(),
                
                fluidRow(
                  
                  br(),
                  column(width=6, h4("Compounds settings"), rHandsontableOutput("Table_compounds_settings")),
                  column(width=4, offset=1, h4("Samples settings"), rHandsontableOutput("Table_samples_settings"))
                  ),
                
                br(),
                  
                fluidRow(column(width=3, offset=9, align="center", 
                                checkboxInput("Checkbox_subtract_MB","Subtract Method Blanks", value=T),
                                shiny::actionButton("Button_generate_heatmap", "Generate Heatmap", icon("chart-bar"), width="200px",
                                                    style="color: #fff; background-color: #2346b0; border-color: #2e6da4"))
                ),
                
                br(), hr(),
                  
                  fluidRow(
                    h3("2. HEATMAP (NORMALIZED)"),
                    br(),
                    checkboxInput("Checkbox_cluster_compounds","Cluster Compounds"),
                    checkboxInput("Checkbox_cluster_samples","Cluster Samples"),
                    br(),
                    h4("log2 fold-change of median-normalized peak areas for each compound"),
                    column(width=12, plotOutput("Plot_heatmap", height="1100px"))
                  ),
                
                br(), br(),
                
                fluidRow(
                  column(width=3, align="center", 
                         fluidRow(
                           shiny::downloadButton("Button_download_normalized_csv", "Normalized Table", icon("file-csv"), width="200px",
                                               style="color: #fff; background-color: #00ab66; border-color: #2e6da4"),
                           shiny::downloadButton("Button_download_normalized_csv_no_qc", "Normalized Table (No QC)", icon("file-csv"), width="200px",
                                               style="color: #fff; background-color: #00ab66; border-color: #2e6da4")
                         )),
                  
                  column(width=3, align="center", shiny::downloadButton("Button_download_heatmap", "Heatmap",
                                                                        icon("file-pdf"), width="200px",
                                                                        style="color: #fff; background-color: #00ab66; border-color: #2e6da4"))
                ),
                
                br(), br()

)

