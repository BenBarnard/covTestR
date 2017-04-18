library(shiny)

shinyUI(fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose File'),
      selectInput('pops', 'Populations', c(2, 3), selectize=FALSE),
      selectInput('xaxis', 'X Axis', c("Differences", "Samples", "Dimensions"), selectize=FALSE),
      uiOutput("plotoptions1"),
      uiOutput("plotoptions2")
    ),
    mainPanel(
      plotOutput(outputId = "plot", height = "300px")
    )
  )
)
)
