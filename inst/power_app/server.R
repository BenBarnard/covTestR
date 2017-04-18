library(shiny)
library(ggplot2)
library(dplyr)

shinyServer(function(input, output, session){

  dataInput <- reactive({
    inFile <- input$file1
    if (is.null(inFile)){
      return(NULL)
    }else{
      file.rename(inFile$datapath,
                  paste0(inFile$datapath, ".RData"))
      load(paste0(inFile$datapath, ".RData"))
      if(is.factor(power_values$Differences)){
        power_values <- mutate(power_values, Differences = as.numeric(as.character(Differences)))
      }
      power_values
    }
  })


   tests <- reactive({
    if (is.null(dataInput())){
      return(NULL)
    }else{
      df <- unique(select(dataInput(), Test))
      df
    }
  })

   differences <- reactive({
     if (is.null(tests())){
       return(NULL)
     }else{
       df <- unique(select(dataInput(), Differences))
       df
     }
   })

   populations <- reactive({
     if (is.null(differences())){
       return(NULL)
     }else{
       df <- unique(select(dataInput(), Populations))
       df
     }
   })

   samples <- reactive({
     if (is.null(populations())){
       return(NULL)
     }else{
       df <- unique(select(dataInput(), Samples))
       df
     }
   })

   dimensions <- reactive({
     if (is.null(samples())){
       return(NULL)
     }else{
       df <- unique(select(dataInput(), Dimensions))
       df
     }
   })



   output$plotoptions1 <- renderUI({
     if(!(is.null(dimensions()))){
       switch(input$xaxis,
              "Differences" = sliderInput("dynamic1", "Samples", min = min(samples()),
                                         max = max(samples()), value = max(samples()),
                                         step = 5),
              "Samples" = sliderInput("dynamic1", "Dimensions", min = min(dimensions()),
                                      max = max(dimensions()), value = max(dimensions()),
                                      step = 5),
              "Dimensions" = sliderInput("dynamic1", "Differences", min = min(differences()),
                                         max = max(differences()), value = max(differences()),
                                         step = 0.1))
     }
   })

   output$plotoptions2 <- renderUI({
     if(!(is.null(dimensions()))){
       switch(input$xaxis,
              "Differences" = sliderInput("dynamic2", "Dimensions", min = min(dimensions()),
                                         max = max(dimensions()), value = max(dimensions()),
                                         step = 5),
              "Samples" = sliderInput("dynamic2", "Differences", min = min(differences()),
                                      max = max(differences()), value = max(differences()),
                                      step = 0.1),
              "Dimensions" = sliderInput("dynamic2", "Samples", min = min(samples()),
                                         max = max(samples()), value = max(samples()),
                                         step = 5))
     }
   })


   output$plot <- renderPlot({
     if(is.null(input$dynamic2)){
       return()
     }else{
       if(input$xaxis == "Differences"){
         df <- select(filter(dataInput(),
                             Samples == input$dynamic1,
                             Dimensions == input$dynamic2,
                             Populations == input$pops),
                      Power, Test, Differences)
         plot <- ggplot(data = df) +
           geom_line(aes(x = Differences, y = Power, color = Test))
       }
       if(input$xaxis == "Samples"){
         df <- select(filter(dataInput(),
                             Dimensions == input$dynamic1,
                             Differences == input$dynamic2,
                             Populations == input$pops),
                      Power, Test, Samples)
         plot <- ggplot(data = df) +
           geom_line(aes(x = Samples, y = Power, color = Test))
       }
       if(input$xaxis == "Dimensions"){
         df <- select(filter(dataInput(),
                             Differences == input$dynamic1,
                             Samples == input$dynamic2,
                             Populations == input$pops),
                      Power, Test, Dimensions)
         plot <- ggplot(data = df) +
           geom_line(aes(x = Dimensions, y = Power, color = Test))
       }
     }
     plot
   })



})
