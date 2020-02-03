tabBox(
  width = 12,
  tabPanel("PCA",
           sidebarLayout(
             sidebarPanel(width = 3,
                          selectInput("y_axisPC", label = "PC on y-axis", ## this should be updated as we figure out how many PCs there are...should be in server, look up how to do this
                                      choices = c('2' = '2'),
                                      selected = '2'),
                          selectInput("x_axisPC", label = "PC on x-axis", ## this should be updated as we figure out how many PCs there are...should be in server, look up how to do this
                                      choices = c('1' = '1'),
                                      selected = "1"),
                          actionButton('genplotpca', 'Generate plot or update plot', icon("chart-bar"),
                                       style="color: #fff; background-color: #006E90; border-color: #006E90"), #337ab7
                          tags$hr(),
                          
                          uiOutput("colourpickers"),
                          tags$hr(),
                          downloadButton('dlPCA', 'Download PCA biplot',
                                         style="color: #fff; background-color: #006E90; border-color: #006E90")
             ),
             mainPanel(
             plotlyOutput("pcaPlot") %>% withSpinner()
             )
)
  ),
  
  tabPanel("Sample clustering",
           sidebarLayout(
             sidebarPanel(width = 3,
                          uiOutput("colourpickers2"),
                          selectInput("dist_method", label = "Distance method:",
                                      choices = c('Euclidean' = 'eucl',
                                                  'Canberra' = 'canb',
                                                  'Binary' = 'bina',
                                                  'Minkowski' = 'mink'),
                                      selected = 'eucl'),
                          selectInput("hclust_method", label = "Hierarchical clustering method:",
                                      choices = c('Ward' = 'ward.D',
                                                  'Ward 2' = 'ward.D2',
                                                  'Single' = 'sing',
                                                  'Complete' = 'compl',
                                                  'Average' = 'aver',
                                                  'McQuitty' = 'mcq',
                                                  'Median' = 'medi',
                                                  'Centroid' = 'centr'),
                                      selected = 'ward.D2'),
                          actionButton('genclustdendro', 'Generate cluster dendrogram',
                                       icon("chart-bar"),
                                       style="color: #fff; background-color: #006E90; border-color: #006E90"),
                          tags$hr(),
                          downloadButton('dlDendro', 'Download cluster dendrogram',
                                         style="color: #fff; background-color: #006E90; border-color: #006E90")
                          ),
             
             mainPanel(
              plotOutput("clustDendro") 
             #plotOutput("clustDendro", height = 900, width = 800)
               )
           )
           
  ),      
  
  tabPanel("Functional enrichment heatmap",
           sidebarLayout(
             sidebarPanel(width = 3,
                          radioButtons("restrict_analysis", "Restrict analysis to compare two conditions?",
                                       c("Yes" = "y", 
                                         "No" = "n"),
                                       selected = "n"),
                          uiOutput("control_gsva"),
                          uiOutput("treatment_gsva"),
                          tags$hr(),
                          radioButtons("fig_type", "Heatmap type:",
                                       c('Classic heatmap' = 'heatmap',
                                         'Bubble heatmap' = 'bubble'),
                                       selected = 'heatmap'),
                          colourInput("high_col", "Colour for high GSVA score", "FF6F59"),
                          colourInput("low_col", "Colour for low GSVA score", "#67A7C1"),
                          tags$hr(),
                          helpText("Please press 'Generate/update heatmap' when you are ready to plot the data. If you would like to update the plot with any additional options or colours, please press the button again."),
                          actionButton('genplotheat', 'Generate/update heatmap',
                                       icon("chart-bar"),
                                       style="color: #fff; background-color: #006E90; border-color: #006E90"),
                          tags$hr(),
                          radioButtons("sample_ord", "Order Samples:",
                                       c('By clustering' = 'clust',
                                         'By conditions/facets' = 'facets')),
                          radioButtons("kegg_ord", "Order KEGG pathways:",
                                       c('By clustering' = 'clust',
                                         'By p-value' = 'pval')),
                          radioButtons('plotsig', "Only plot significantly enriched KEGG?",
                                       c('Yes' = 'y',
                                         'No' = 'n'),
                                       selected = 'y'),
                          conditionalPanel(
                            condition = "input.plotsig == 'y'",
                            textInput("pvalthresh", "Adjusted p-value threshold for plotting", "0.05"))
#                          )
,

                          
                          tags$hr(),
                          downloadButton('downloadPlot','Download heatmap',
                                         style="color: #fff; background-color: #006E90; border-color: #006E90"),
                          tags$hr(),
                          downloadButton('downloadKEGG', 'Download peptide annotation',
                                         style="color: #fff; background-color: #006E90; border-color: #006E90")
             ),
             mainPanel(
               #plotOutput("heatmapPlot")
              
               
               fluidRow(
                 tabBox(
                   title = tagList(shiny::icon("gear"), "Plot dimension settings"), width = 12,
                   tabPanel("In app visualization", 
                            sliderInput("plotheight", "Plot height (pixels):",
                             min = 100, max = 1000, value = 500),
                            sliderInput("plotwidth", "Plot width (pixels):",
                            min = 100, max = 1000, value = 700)
                    ),
                   tabPanel("To PDF",
                            sliderInput("plotheightsave", "Plot height (inches):",
                                         min = 1, max = 20, value = 8),
                            sliderInput("plotwidthsave", "Plot width (inches):",
                                      min = 1, max = 20, value = 11)
                            )
                 )
                 
               ),
               fluidRow(
                 box(width=12,
                     uiOutput("heatmapUI")
                 )
               ),

           )
  )
))
