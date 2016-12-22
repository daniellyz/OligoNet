library(shiny)
library(KEGGREST)
require(visNetwork, quietly = TRUE)
library(markdown)

textInputRow<-function (inputId, label, value = "") 
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}
kegg_organism=read.table("kegg_organisms.txt",sep="\t",header=F,stringsAsFactors = F)


shinyUI(fluidPage(
  
  titlePanel('OligoNet: Oligopeptides in Metabolomics Data'),
  
  # Load D3.js
  tags$head(
    tags$script(src = 'http://d3js.org/d3.v3.min.js')
  ),
  
  fluidRow(
    column(4,
         h3('A) Upload metabolomics data'),
         p('Example of file 1:', a(href = 'https://raw.githubusercontent.com/daniellyz/Curl/master/matrix.txt', 'matrix.txt')),
         p('Example of file 2:', a(href = 'https://raw.githubusercontent.com/daniellyz/Curl/master/matrix_additional.txt', 'matrix_additional.txt')),
         
         fileInput('file1', 'File 1: ID + Neutral Mass + Intensities',
                   accept = c('dat/tab-separated-value','text/tab-separated-values')),
         fileInput('file2', 'File 2: ID + Additional Informations (Optional)',
                   accept = c('dat/tab-separated-value','text/tab-separated-values'))),

    column(4,
        h3('B) Perform decomposition'),
           
        p('Default alphabets (subunits):', a(href = 'https://raw.githubusercontent.com/daniellyz/Curl/master/amino-acid.txt', 'amino-acid.txt')),
        fileInput('file3', label=h6('File 3: modify the default amino acid file (optional)'),
        accept = c('dat/tab-separated-value','text/tab-separated-values')),
        radioButtons("TE", label="Please select the nature of your mass:", 
        choices = list("Theoritical mass" = 1, "Experimental mass with error (Da):" = 2),selected = 1),
        numericInput('tol', h6(''), 0.01,min=1e-5,max=0.5),
        textInput("DecompID", label = "Enter Bibiserv2-Decomp job ID (Optional):",value="Job ID"),
        actionButton("goButton", "Go!"),
        downloadButton("downloadAnnotation", "Download decomposition results")),   
         
    column(4,
        h3('C) Network construction'),
        checkboxGroupInput("visual", label = h6(""), choices = list("Show entire network" = 1, "Remove free amino acids" = 2, "Correlation between nodes bigger than" = 3), selected = 1),
        sliderInput("cor_min", label = h6(""), min = -1, max = 1, value = -1,step=0.1),
        uiOutput("Shows"),
        actionButton("goButtonbis", "Go!"),
        downloadButton("downloadNetwork", "Download network (edges) file"))),
    
  
  br(),
  br(),
  
    # Show network graph
    mainPanel(
      tabsetPanel(
        tabPanel("Help",includeMarkdown("Help.Rmd")),
        tabPanel("Job Status",verbatimTextOutput("summary")),
        tabPanel("UAAC-Annotation",dataTableOutput("table1")),
        tabPanel("MAAP-Annotation",dataTableOutput("table2")),
        tabPanel("Entire PDN",visNetworkOutput("networkPlot",height = "1000px")),
        tabPanel("PDN analysis",fluidRow(
                                  h3('Degree distribution:'),
                                  plotOutput("Distribution_degree",height =300, width = 550),
                                  h3('Edge distribution:'),
                                  plotOutput("Distribution_edge",height =300, width = 550),
                                  h3('Top 20% most frequent edges:'),
                                  plotOutput("Top_edge",height =300, width = 2000),
                                  h3('Path length distribution:'),
                                  plotOutput("Distribution_length",height =300, width = 550),
                                  h3('Spearman correlation distribution:'),
                                  plotOutput("Distribution_correlation",height =300, width = 550))),
        tabPanel("Subgraphs", fluidRow(
                                  br(),
                                  h3('High degree vertices:'),
                                  uiOutput("Controls"),
                                  selectInput('Order', 'Neighbourhood order:',1:5),
                                  visNetworkOutput("networkPlot3",height = "500px"),
                                  downloadButton("downloadNetwork3", "Download network (edges) file"),
                                  h3('Degradation chains:'),
                                  br(),
                                  uiOutput("Chains"),
                                  visNetworkOutput("networkPlot1",height = "500px"),
                                  downloadButton("downloadNetwork1", "Download network (edges) file"))),
        tabPanel("Peptides in metabolic pathways ",fluidRow(
                                br(),
                                h3('ATTENTION! This function can be time-consuming'),
                                br(),
                                column(5,
                                selectInput('Organism', 'Choose your organism:',kegg_organism[,2])),
                                column(5,offset = -5,align="bottom",
                                uiOutput("Kegg_pathway")),
                                br(),
                                actionButton("goButton3", "Go!")
                                ),
                 htmlOutput('image')),
        tabPanel("About",includeMarkdown("About.Rmd"))
       )
    )

))


