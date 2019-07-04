options(encoding = "UTF-8")
library("shiny")
library("ggplot2")
library("phyloseq")

path = "D:/stage/data/runs_new2"
setwd(path)
load("objects.Rdata")



shinyUI( pageWithSidebar( 
  
  headerPanel("Test"), # titre de l'appli
  
  sidebarPanel( # cette partie va contenir tous les éléments de contrôle de l'UI
    
    selectInput('dataset', 'Choose a dataset:', choices = c("ps", 
                                                            "ps_decontam", 
                                                            "ps_decontam2",
                                                            "ps_percent",
                                                            "ps_proteo",
                                                            "ps_wolbachia")), 
    
    radioButtons('format', "Choose a output format:", choices =c("png",
                                                                 "jpeg", 
                                                                 "pdf")),
    
    downloadButton("downloadPlot")
  ),
  

  
  mainPanel( # cette partie va contenir les sorties
    
    h3("Résultats :"), # titre donné à la partie présentant les sorties, l'élément h3 correspond à une balise "<h3>" en html, ie. le titre sera mis en valeur
    verbatimTextOutput("summary"),
    plotOutput("mon_plot") 
    
  )
  
))
