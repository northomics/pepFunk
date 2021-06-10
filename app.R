library(rhandsontable)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(colourpicker)
library(reshape2)
library(DT)
library(tidyverse)
library(plyr)
library(DESeq2)
library(GSVA)
library(limma)
library(ggdendro)
library(plotly)
library(dendextend)
library(LaCroixColoR) #devtools::install_github("johannesbjork/LaCroixColoR")
library(shinycssloaders)


# install.packages.auto(rhandsontable)
# devtools::install_github("johannesbjork/LaCroixColoR")
# install.packages.auto(DESeq2)
# install.packages.auto(GSVA)
#
# # if (!requireNamespace("BiocManager", quietly = TRUE))
# # install.packages("BiocManager")
# #
# # BiocManager::install("DESeq2")
# # BiocManager::install("GSVA")

# The proteinGroup files are large, increase the limit of file upload size to 500MB
options(shiny.maxRequestSize=500*1024^2)

source("peptide_centric_functions.R")
source("peptide_centric_module.R")



########
# DATA # # should add files to this folder
########

## Full KEGG pathways database
kegg_L3 <- read.delim("./data/kegg_L3.tsv", sep='\t', header=F,
                      col.names=c('L3', 'L3_desc', 'L4', 'L4_desc'),
                      colClasses=c('character','character','character','character')) %>% as.data.frame()
KOpathways <- kegg_L3$L4_desc %>% unique()
pathway_kegg <- dlply(kegg_L3 %>% dplyr::select(L4_desc, L3), .(L4_desc))

## Full COG pathways database
cog_categories <- read.delim("./data/COG.2020.categories.tsv", sep='\t', header=F,
                      col.names=c('L3', 'L3_desc', 'L4', 'L4_desc','L5'),
                      colClasses=c('character','character','character','character','character')) %>% as.data.frame()
COGpathways <- cog_categories$L4_desc %>% unique()
pathway_COG <- dlply(kegg_L3 %>% dplyr::select(L4_desc, L3), .(L4_desc))

## peptide to KEGG database now loaded depending on user's choice in server_analyze.R

#  _header ------------------------------------------------------


header <- dashboardHeader(title = span(img(src="logo_imetalab.png", width = 140), "pepFunk"),
                          titleWidth = 460,
                          tags$li(class = "dropdown",
                                  tags$a(tags$img(height = "18px", alt="SNAP Logo", src="logo_M.png")
                                  ),
                            tags$head(
                                    tags$link(rel = "stylesheet", type = "text/css", href = "./www/custom.css")
                                  )
                          )
)


#  _side bar ------------------------------------------------------

sidebar <- dashboardSidebar(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),
  width = 250,
  sidebarMenu(
    id = "tabs",
    menuItem("Upload Data", tabName = "dashboard", icon = icon("table")),
    # Get icon codes from here: https://fontawesome.com/v4.7.0/icons/
    sidebarMenuOutput("menu1"),
    menuItem("Analysis", tabName = "Analysis", icon = icon("cog")),
    menuItem("Gallery", tabName = "gallery", icon = icon("picture-o")),
    menuItem("About", tabName = "about", icon = icon("question-circle")),
    menuItem("iMetaLab", icon = icon("home"),
             href = "http://www.imetalab.ca"),
    menuItem("pepFunk on GitHub", icon=icon("github"),
             href = "https://github.com/northomics/pepFunk")

  )
)


# _body --------------------------------------------------------------

body <- dashboardBody(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),
  tabItems(
    #  ___dashboard/starting  tab  ------------------------------------------------------
    tabItem(tabName = "dashboard",
            source(
              file = "ui_dashboard.R",
              local = TRUE,
              encoding = "UTF-8"
            )$value),
    tabItem(tabName = "Analysis",
            source(
              file = "ui_analyze.R",
              local = TRUE,
              encoding = "UTF-8"
            )$value),
    tabItem(tabName = "gallery",
            tabItem(tabName = "gallery",
                    fluidRow(
                      box(
                        title = "PCA analysis",
                        solidHeader = TRUE,
                        status = "primary",
                        width = 6,
                        img(src='gallery/PCA.png')
                      ),
                      box(
                        title = "Heatmap of enriched functions",
                        solidHeader = TRUE,
                        status = "primary",
                        width = 6,
                        img(src='gallery/Heat.png')
                      )
                    )
            )
            ),
    tabItem(tabName = "about",
            fluidRow(
              box(
                title = "About pepFunk, a metaproteomic peptide-centric functional enrichment workflow.",
                solidHeader = TRUE,
                status = "primary",
                width = 12,
                "Welcome to pepFunk!",
                br(),
                "pepFunk allows you to complete a peptide-focused functional enrichment workflow for gut microbiome metaproteomic studies.

                This workflow uses KEGG, COG, or eggNOG annotation for pathway enrichment, alongside Gene Set Variation Analysis (GSVA) adapted for peptide data.
                By completing analysis on peptides, rather than proteins, we lose less information and retain more statistical power.
                We curated peptide database specific to human gut microbiome studies for computational speed.
              "
              )
            ),
            fluidRow(
              box(
                title = "Updates",
                solidHeader = TRUE,
                status = "primary",
                width = 12,
                "January 22, 2021",
                br(),
                "Fixed a bug where manual condition entering would cause an error.
              ",
                br(),
                "February 22, 2021",
                br(),
                "Fixed a bug to allow 10 conditions."

              )
            )
    )
  ),
  #Semi-collapsible sidebar
  tags$script(HTML("$('body').addClass('sidebar-mini');"))
)



# ------ UI ---------------------------
ui <- dashboardPage(
  title = "pepFunk",
  header,
  sidebar,
  body,
)


server <- function(input, output, session) {
  source(file = "server_analyze.R",
         local = TRUE,
         encoding = "UTF-8")
  output.fileUploaded1 = TRUE
}

shinyApp(ui, server)
