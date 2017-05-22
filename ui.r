library(shiny)
library(KEGGREST)
require(visNetwork, quietly = TRUE)
library(markdown)
library(DT)

textInputRow<-function (inputId, label, value = "") 
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}
kegg_organism=read.table("kegg_organisms.txt",sep="\t",header=F,stringsAsFactors = F)


shinyUI(navbarPage("OligoNet: Oligopeptides in Metabolomics Data",

  # Load D3.js
 # tags$head(
#   tags$script(src = 'http://d3js.org/d3.v3.min.js'),
   # tags$style("label {display:inline;}")
    
  #),
  
 # br(),
  #br(),
  
  # Show network graph

    tabPanel("A) Start a run",
     column(6,         
        br(),
        em('* means the information is compulsory'),
        br(),
        h3("Upload metabolomics data (IDs + Mass signals+ Intensities)*"),
        fileInput('file1',label='',accept = c('.txt','.dat','.csv')),
        h4("or paste your data into the field below:"), 
        textAreaInput("blank_file1", label = '',width=500,height=200),
        br(),
        h3("Upload metabolomics data (IDs + Other information)"),
        fileInput('file2',label='',accept = c('.txt','.dat','.csv')),
        #https://github.com/daniellyz/OligoNet/blob/master/amino-acid-basic.txt'
        h3('Modify the default monomer (Amino acid) file'),
        fileInput('file3', label='',accept = c('.txt','.dat','.csv')),
        h3("Scan mode*:"),
        selectInput('Scan', label='',c('Neutral (input data is already corrected)','Positive ionization (correct for H+ adduct)','Negative ionization (correct for H+ loss)')),
        h3("Max error (Da, 0 Da for theoritical masses)*:"),
        numericInput('tol', label='',0,min=0,max=0.5,step=0.001),
        h3("Enter Bibiserv2-Decomp job ID:"),
        textInput("DecompID", label='',value="Job ID"),
        h3("Faster peptide annotation (For more than 1000 mass features):"),
        br(),
        checkboxInput("checkbox", label = "Use DECOMP server (faster)", value = TRUE),
        br(),
        br(),
        actionButton("goButton", "Start the run",style='padding:6px; font-size:150%')
     ),
    
    column(5,
        em('Messages from the server:'),
        verbatimTextOutput("blank_message1")
    )),

    tabPanel("B) Mass decomposition",
        
        br(),
        downloadButton("downloadAnnotation", "Download decomposition results",style='padding:6px; font-size:150%'),
        br(),
        br(),
        
        tabsetPanel("",
        tabPanel("UAAC-Annotation",
          h3("Mass signals annotated to Unique Amino Acid Combination"),
          br(),
          dataTableOutput("table1")),
        
        tabPanel("MAAP-Annotation",
          h3("Mass signals annotated to Multiple Amino Acid Combination"),
          br(),
          dataTableOutput("table2")))),
        
        
    tabPanel("C) Network construction",
        
        br(),
        h3("Please choose one or more options*:"),
        
        HTML("<div style='height: 60px;'>"),
        checkboxGroupInput("visual", label = h6(""), choices = list("Show entire network" = 1, "Remove triangles" = 2,"Remove free amino acids" = 3), selected = c(1,2)),
        HTML("</div>"),
        div(style="display: inline-block;vertical-align:middle; width: 300px;",checkboxGroupInput("visual2", label = h6(""), choices =list("Correlation between nodes" = 4), selected = 0)),
        div(style="display: inline-block;vertical-align:middle; width: 300px;",sliderInput("cor_min_max", label = h6(""), min = -1, max = 1, value = c(0.8,1),step=0.1)),
        uiOutput("Shows"),
        uiOutput("Coloring"),
        uiOutput("Coloring_nodes"),
        
        br(),
        actionButton("goButtonbis", "Start network construction",style='padding:6px; font-size:150%')),


    tabPanel("Network results",        
        
        br(),     
        downloadButton("downloadNetwork", "Download network (edges) file",style='padding:6px; font-size:150%'),      
        br(),
        br(),
        
        tabsetPanel("",
         tabPanel("Entire network",
            h3("Peptide Degradation Network"),
            visNetworkOutput("networkPlot",height = "1000px")),
        
         tabPanel("Network analysis",
           h3('Degree distribution:'),
           plotOutput("Distribution_degree",height =300, width = 400),
           h3('Edge distribution:'),
           plotOutput("Distribution_edge",height =300, width = 400),
           h3('Top 20% most frequent edges:'),
           plotOutput("Top_edge",height =300, width = 1200),
           h3('Path length distribution:'),
           plotOutput("Distribution_length",height =300, width = 400),
           h3('Spearman correlation distribution:'),
           plotOutput("Distribution_correlation",height =300, width = 450)),
         
        tabPanel("Subgraphs", 
          h3('High degree vertices:'),
          uiOutput("Controls"),
          selectInput('Order', 'Neighbourhood order:',1:5),
          visNetworkOutput("networkPlot3",height = "500px"),
          downloadButton("downloadNetwork3", "Download network (edges) file",style='padding:6px; font-size:150%'),
          h3('Degradation chains:'),
          uiOutput("Chains"),
          visNetworkOutput("networkPlot1",height = "500px"),
          downloadButton("downloadNetwork1", "Download network (edges) file",style='padding:6px; font-size:150%')
      ))),
        
      tabPanel("Peptides in KEGG",
        br(),
        h3("Peptides in metabolic pathways"),
        selectInput('Organism', 'Choose your organism*:',kegg_organism[,2]),
        uiOutput("Kegg_pathway"),
        htmlOutput('image')),
      
      tabPanel("Example datasets"),
      tabPanel("Help",includeMarkdown("Help.Rmd")),
      tabPanel("About",includeMarkdown("About.Rmd"))
      )
      
  
)


