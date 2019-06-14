options(encoding = "UTF-8")
library("shiny")
library("ggplot2")
library("phyloseq")

load("objects.Rdata")



shinyServer( function(input, output) { # les éléments de la partie 'server' vont être définis ici
  
  
  datasetInput <- reactive({
    switch(input$dataset,
           "ps" = ps,
           "ps_decontam" = ps_decontam,
           "ps_decontam2" = ps_decontam2,
           "ps_percent" = ps_percent,
           "ps_proteo" = ps_proteo,
           "ps_wolbachia" = ps_wolbachia)
  })
  
  output$summary <- renderPrint({
    dataset <- datasetInput()
    dataset
  })
  

  
  
  plot_input <- reactive({ 
    mydata <- datasetInput()
    
    readsumsdf <- data.frame(nreads = sort(taxa_sums(mydata), TRUE),
                             sorted = 1:ntaxa(mydata), 
                             type = "OTU")
    
    
    readsumsdf2 <- data.frame(nreads = sort(sample_sums(mydata), TRUE), 
                              sorted = 1:nsamples(mydata), 
                              type = "Samples")
    
    readsumsdf3 <- rbind(readsumsdf,readsumsdf2)
    
    
    p  <-  ggplot(readsumsdf3, 
                  aes(x = sorted, y = nreads)) + 
      geom_bar(stat = "identity")+
      theme_gray() +
      ggtitle("Total number of reads before Preprocessing") + 
      scale_y_log10() +
      facet_wrap(~type, 1, scales = "free")
    
  }) 

    
    
output$mon_plot <- renderPlot({
    print(plot_input())
  })
  

  output$downloadPlot <- downloadHandler(
    filename = function(){
      paste("distribution",input$format, sep='.')
    },
    content = function(file) {
      if(input$format == "jpeg"){
        ggsave(file, plot = plot_input(), device = "jpeg")
      }
      if (input$format == "pdf"){
        ggsave(filename, plot = plot_input(), device = "pdf")
      }
      if (input$format == "png"){
        ggsave(file, plot = plot_input(), device = "png")
      }
    }
  )
  
})