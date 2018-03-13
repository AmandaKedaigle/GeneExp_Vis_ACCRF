library(shiny)
library(ggplot2)
library(datasets)

{ ##DATA PREP
  load("expression_data.RData")
  tpm = within(tpm, {
    sample = ifelse(variable %in% c("ERR315325.bam","ERR315382.bam","ERR315418.bam","ERR315420.bam","ERR315449.bam","ERR315459.bam"), "Normal Salivary Gland", "ACC tumor")
  })
  #colnames(marray) = c("Transcript.Cluster.ID","ACC Tumor","Normal Salivary Gland","FoldChange","ANOVA.p","FDR","Gene.Symbol","Description")
  #mmarray <- melt(marray[,c(1,2,3,7)],id.vars=c("Gene.Symbol","Transcript.Cluster.ID"))
}

# Define server logic
shinyServer(function(input, output, session) {
  
  updateSelectizeInput(session, 'gene', choices=tpm$Gene, server=TRUE)
  
  output$rnaseqPlot <- renderPlot({
    ggplot(tpm[tpm$Gene==input$gene,], aes(variable,value,color=sample)) +
      geom_point(size=5) +
      xlab(NULL) + ylab("Transcripts per Million") +
      labs(title="RNA-seq") +
      theme(axis.text.x = element_text(angle = 50, hjust = 1))
  })
  
  output$microarrayPlot <- renderPlot({
    if (input$gene %in% mmarray$Gene.Symbol) {
      ggplot(mmarray[mmarray$Gene.Symbol==input$gene,], aes(variable,value,color=variable)) +
        geom_line(aes(group = Transcript.Cluster.ID),color="black") +
        geom_point(size=5) +
        xlab(NULL) + ylab("Bi-weight Avg Signal (log2)") +
        labs(title="Microarray", 
             caption=paste0("FoldChange = ",marray[marray$Gene.Symbol==input$gene,"FoldChange"],", FDR = ",marray[marray$Gene.Symbol==input$gene,"FDR"])) +
        theme(axis.text.x = element_text(angle = 50, hjust = 1),
              plot.caption = element_text(size=14))
    } else {
      plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
      text(x=0.4, y=1, "No microarray data for this gene", cex = 1.6, col = "black")
    }
  })
  
})
