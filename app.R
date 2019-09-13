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
library(dendextend)
library(viridis)
library(LaCroixColoR) #devtools::install_github("johannesbjork/LaCroixColoR")

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

# The proteinGroup files are large, increnase the limit of file upload size to 500MB
options(shiny.maxRequestSize=500*1024^2)

source("peptide_centric_functions.R")
source("peptide_centric_module.R")



########
# DATA # # should add files to this folder
########

## Full KEGG database
kegg_L3 <- read.delim("./www/kegg_L3.txt", sep='\t', header=F, 
                      col.names=c('L3', 'L3_desc', 'L4', 'L4_desc'),
                      colClasses=c('character','character','character','character')) %>% as.data.frame()
pathways <- kegg_L3$L4_desc %>% unique() 
pathways <- pathways[-c(208:229)] #removing KO that are not in brite/pathway
pathway_kegg <- dlply(kegg_L3 %>% dplyr::select(L4_desc, L3), .(L4_desc))

## Core kegg database
core_pep_kegg <- read.delim("./www/core_pep_kegg.csv", 
                            sep=",", header=F, col.names = c("pep", "kegg", "count", "eval"))
core_pep_kegg_only <- core_pep_kegg %>% dplyr::group_by(pep) %>% dplyr::select(pep, kegg)
core_pep_kegg <- core_pep_kegg %>% dplyr::group_by(pep) %>%
  dplyr::summarize(total = sum(count))  %>%
  merge(., core_pep_kegg, by='pep', all.y=T) %>%
  dplyr::mutate(prop = count/total) %>% dplyr::select(pep, kegg, prop)
core_pep_kegg$newpep_name <- make.names(core_pep_kegg$pep,unique=T) # update pep names so all unique

#  _header ------------------------------------------------------


header <- dashboardHeader(title = span(img(src="logo_imetalab.png", width = 140), "Peptide-Centric Analyst"),
                          titleWidth = 460,
                          tags$li(class = "dropdown",
                                  tags$a(tags$img(height = "18px", alt="SNAP Logo", src="logo_M.png")
                                  )
                          )
)


#  _side bar ------------------------------------------------------

sidebar <- dashboardSidebar(
  width = 250,
  sidebarMenu(
    id = "tabs",
    menuItem("Upload Data", tabName = "dashboard", icon = icon("file")),
    # Get icon codes from here: https://fontawesome.com/v4.7.0/icons/
    sidebarMenuOutput("menu1"),
    menuItem("Analysis", tabName = "Analysis", icon = icon("file")),
    menuItem("Gallery", tabName = "gallery", icon = icon("picture-o")),
    menuItem("About", tabName = "about", icon = icon("question-circle")),
    menuItem("iMetaLab", icon = icon("home"), 
             href = "http://www.imetalab.ca")
    
  )
)


# _body --------------------------------------------------------------

body <- dashboardBody(
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
                title = "About protein-centric metaproteomic workflows...", 
                solidHeader = TRUE,
                status = "primary", 
                width = 12,
                "Welcome to the Peptide-Centric Analyst!",
                br(),
                "Here you can write down information about your app"
              )
            )   
    )
  ),

  #   CSS section ignored for analysis------------------------------------------------------ 
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),
  #Semi-collapsible sidebar
  tags$script(HTML("$('body').addClass('sidebar-mini');"))             
)



# ------ UI ---------------------------
ui <- dashboardPage(
  title = "Peptide-Centric Analyst",
  header,
  sidebar,
  body
)


server <- function(input, output, session) {
  source(file = "server_analyze.R",
         local = TRUE,
         encoding = "UTF-8")
  output.fileUploaded1 = TRUE
}

shinyApp(ui, server)