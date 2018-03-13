library(shiny)

# Use a fluid Bootstrap layout
shinyUI(fluidPage(    
  
  # Give the page a title
  titlePanel("Gene Expression in ACC PDX models"),
  
  # Generate a row with a sidebar
  sidebarLayout(      
    
    # Define the sidebar with one input
    sidebarPanel(
      selectizeInput("gene", "Gene:", choices=NULL,options=list(placeholder="e.g., 'NOTCH1'")),
      hr(),
      helpText("Enter a gene symbol to view mRNA expression levels")
    ),
    
    # Create a spot for the scatter plot
    mainPanel(
      plotOutput("rnaseqPlot"), 
      plotOutput("microarrayPlot")
    )
    
  )
)
)
