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
                          actionButton('genplotpca', 'Generate plot and update colours'),
                          tags$hr(),
                          
                          uiOutput("colourpickers"),
                          tags$hr(),
                          downloadButton('dlPCA', 'Download PCA biplot')
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
                          actionButton('genclustdendro', 'Generate cluster dendrogram'),
                          tags$hr(),
                          downloadButton('dlDendro', 'Download cluster dendrogram')
                          ),
             
             mainPanel(plotlyOutput("clustDendro"))
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
                          colourInput("high_col", "Colour for high GSVA score", "FF6F59"),
                          colourInput("low_col", "Colour for low GSVA score", "#67A7C1"),
                          radioButtons("sample_ord", "Order Samples:",
                                       c('By clustering' = 'clust',
                                         'By conditions/facets' = 'facets')),
                          radioButtons("kegg_ord", "Order KEGG pathways:",
                                       c('By clustering' = 'clust',
                                         'By p-value' = 'pval')),
                          radioButtons('plotsig', "Only plot significantly enriched KEGG?",
                                       c('Yes' = 'y',
                                         'No' = 'n'),
                                       selected = 'n'),
                          conditionalPanel(
                            condition = "input.plotsig == 'y'",
                            textInput("pvalthresh", "Adjusted p-value threshold for plotting", "0.05")),
                          actionButton('genplotheat', 'Generate/update heatmap'),
                          downloadButton('downloadPlot','Download heatmap'),
                          tags$hr(),
                          downloadButton('downloadKEGG', 'Download peptide annotation')
             ),
             mainPanel(
               plotlyOutput("heatmapPlot")
             )
           )
  )
)
