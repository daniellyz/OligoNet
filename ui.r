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
         
         fileInput('file1', 'File 1: ID + Neutral Mass + Intensities',
                   accept = c('dat/tab-separated-value','text/tab-separated-values')),
         fileInput('file2', 'File 2: ID + Additional Informations (Optional)',
                   accept = c('dat/tab-separated-value','text/tab-separated-values')),
         
         actionButton("Example1", "Run Example 1: Yeast metabolic profiling by FT-ICR-MS"),
         h3(),
         actionButton("Example2", "Run Example 2: Yeast metabolic profiling by UPLC-MS"),
         h3(),
         actionButton("clearButton", "Clear examples")
         ),

    column(4, offset = 0,
        h3('B) Perform decomposition'),
           
        checkboxInput("checkbox", label = "Use DECOMP server (Faster)", value = FALSE),
       # p('Default alphabets (subunits):', a(href = 'https://github.com/daniellyz/OligoNet/blob/master/amino-acid-basic.txt', 'amino-acid.txt')),
        fileInput('file3', label=a(href = 'https://github.com/daniellyz/OligoNet/blob/master/amino-acid-basic.txt','File 3: Modify the default amino acid file (optional)'),
        accept = c('dat/tab-separated-value','text/tab-separated-values')),
        radioButtons("TE", label="Please select the nature of your mass:", 
        choices = list("Theoritical mass" = 1, "Experimental mass with error (Da):" = 2),selected = 1),
        numericInput('tol', '', 0.01,min=1e-5,max=0.5),
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
        tabPanel("UAAC-Annotation",fluidRow(
                     h3("Mass signals annotated to Unique Amino Acid Combination"),
                     h4("Help: The column \"Peptide\" contains the decomposition results of each mass signal. Rule: a mass signal that is decomposed into two alanines and two prolines is annotated as A2P2.
                       The column \"KEGG\" represent the KEGG compound code(s) if the mass signal is also annotated in the KEGG databse. Clicking on the link will take you to
                        the webpage describing the compound(s)." ),
                     br(),
                     dataTableOutput("table1"))),
        tabPanel("MAAP-Annotation", fluidRow(
                 h3("Mass signals annotated to Multiple Amino Acid Combination"),
                 h4("Help: The column \"Peptide\" contains the decomposition results of each mass signal. The possible combinations are separated by a semicolon. The column \"NBP\" indicates the number of possible combinations.
                     The column \"KEGG\" represent the KEGG compound code(s) if the mass signal is also annotated in the KEGG databse. Clicking on the link will take you to
                        the webpage describing the compound(s)." ),
                 dataTableOutput("table2"))),
        tabPanel("Entire PDN",fluidRow(h3("Peptide Degradation Network"),
                                       h4("The entire Peptide Degradation Network built from UAAC signals can be visualized:"),         
                                       visNetworkOutput("networkPlot",height = "1000px"))),
        tabPanel("PDN analysis",fluidRow(
                                  h3('Degree distribution:'),
                                  h4("Help: The degree distribution (in + out degree) of all vertices in the PDN"),
                                  plotOutput("Distribution_degree",height =300, width = 550),
                                  h3('Edge distribution:'),
                                  h4("Help: The frequency of edges (different types of degradation reactions) in the PDN"),
                                  plotOutput("Distribution_edge",height =300, width = 550),
                                  h3('Top 20% most frequent edges:'),
                                  plotOutput("Top_edge",height =300, width = 2000),
                                  h3('Path length distribution:'),
                                  h4("Help: All paths in the PDN are extracted. Here we plot the distribution of their length. For instance, the path
                                     \"A2P2 %->% A1P2 %->% P2%->% P1\" has a length 3."),
                                  plotOutput("Distribution_length",height =300, width = 550),
                                  h3('Spearman correlation distribution:'),
                                  h4("Help: The Spearman correlation coefficients between connected vertices in the PDN are calculated. Here is their distribution:"),
                                  plotOutput("Distribution_correlation",height =300, width = 550))),
        tabPanel("Subgraphs", fluidRow(
                                  h3('High degree vertices:'),
                                  h4("Help: Here we list all vertices that have an in-degree of at least 5. Users could choose which vertice they want to visualize:"),
                                  uiOutput("Controls"),
                                  h4("Help: Users could also choose the neighbouring nodes of order 1 to 5 of the selected vertice"),
                                  selectInput('Order', 'Neighbourhood order:',1:5),
                                  visNetworkOutput("networkPlot3",height = "500px"),
                                  downloadButton("downloadNetwork3", "Download network (edges) file"),
                                  br(),
                                  h3('Degradation chains:'),
                                  h4("Help: Here we list all paths longer than 3 in the PDN. Users could choose the path they want to visualize according to the starting node of the path."),
                                  uiOutput("Chains"),
                                  visNetworkOutput("networkPlot1",height = "500px"),
                                  downloadButton("downloadNetwork1", "Download network (edges) file"))),
        tabPanel("Peptides in metabolic pathways ",fluidRow(
                                h3("Peptides in metabolic pathways"),
                                h4("Help:KEGG pathway maps are generated according to the selected organism. In the selected pathway map, annotated peptides (both UAAC and MAAP signals) are annotated in Red, 
                            while other annotated metabolites in the metabolomic dataset are colored in Yellow. Non-annoated metabolites are in white."),
                                selectInput('Organism', 'Choose your organism:',kegg_organism[,2]),
                                uiOutput("Kegg_pathway"),
                                htmlOutput('image'))),
        tabPanel("About",includeMarkdown("About.Rmd"))
       )
    )

))


