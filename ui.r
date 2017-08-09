library(shiny)
library(KEGGREST)
require(visNetwork, quietly = TRUE)
library(markdown)
library(DT)
library("V8")
library(shinyjs)

textInputRow<-function (inputId, label, value = "") 
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}
kegg_organism=read.table("kegg_organisms.txt",sep="\t",header=F,stringsAsFactors = F)


shinyUI(navbarPage("OligoNet: Oligopeptides in Metabolomics Data",

  # Load D3.js
  #tags$head(
  # tags$script(src = 'http://d3js.org/d3.v3.min.js'),
  #  tags$style("label {display:inline;}")
  #),
  
# Show network graph

    tabPanel("A) Start a run",
     shinyjs::useShinyjs(),
     shinyjs::extendShinyjs(text = "shinyjs.refresh = function() { location.reload(); }"),
             
     column(7,         
        br(),
        h3("Upload Compulsory File 1: IDs + Mass signals+ Intensities"),
        div(style="display: inline-block;vertical-align:top; width: 200px;",h5("Please check the examples:")),
        div(style="display: inline-block;vertical-align:top; width: 40px;",h5(a("Ex 1,",href = 'https://raw.githubusercontent.com/daniellyz/OligoNet/master/Datasets/Yeast-FT-Pos.txt',target="_blank"))),
        div(style="display: inline-block;vertical-align:top; width: 40px;",h5(a("Ex 2,",href = 'https://raw.githubusercontent.com/daniellyz/OligoNet/master/Datasets/Yeast-LC-Pos.txt',target="_blank"))),
        div(style="display: inline-block;vertical-align:top; width: 40px;",h5(a("Ex 3,",href = 'https://raw.githubusercontent.com/daniellyz/OligoNet/master/Datasets/Wine_data.txt',target="_blank"))),
        div(style="display: inline-block;vertical-align:top; width: 40px;",h5(a("Ex 4",href = 'https://raw.githubusercontent.com/daniellyz/OligoNet/master/Datasets/Positive_Doping.txt',target="_blank"))),
        
        fileInput('file1',label='',accept = c('.txt','.dat','.csv')),
        h4("or paste your data into the field below:"), 
        textAreaInput("blank_file1", label = '',width=500,height=200),
        br(),
        
        h3("Upload Optional File 2: IDs + Other informations"),
        div(style="display: inline-block;vertical-align:top; width: 200px;",h5("Please check the examples:")),
        div(style="display: inline-block;vertical-align:top; width: 40px;",h5(a("Ex 1,",href = 'https://raw.githubusercontent.com/daniellyz/OligoNet/master/Datasets/Yeast-FT-Pos_additional.txt',target="_blank"))),
        div(style="display: inline-block;vertical-align:top; width: 40px;",h5(a("Ex 2,",href = 'https://raw.githubusercontent.com/daniellyz/OligoNet/master/Datasets/Yeast-LC-Pos-additional.txt',target="_blank"))),
        div(style="display: inline-block;vertical-align:top; width: 40px;",h5(a("Ex 3",href = 'https://raw.githubusercontent.com/daniellyz/OligoNet/master/Datasets/Wine_data_add.txt',target="_blank"))),
        fileInput('file2',label='',accept = c('.txt','.dat','.csv')),
        #https://github.com/daniellyz/OligoNet/blob/master/amino-acid-basic.txt'
       
        h3('Modify the default monomer (Optional File 3: Amino acids)'),
        div(style="display: inline-block;vertical-align:top; width: 260px;",h5("Please check the default amino acid file:")),
        div(style="display: inline-block;vertical-align:top; width: 80px;",h5(a("Amino acid",href = 'https://raw.githubusercontent.com/daniellyz/OligoNet/master/amino-acid-complete.txt',target="_blank"))),
        fileInput('file3', label='',accept = c('.txt','.dat','.csv')),
        h3("Scan mode:"),
        selectInput('Scan', label='',c('Neutral (input data is already corrected)','Positive ionization (correct for H+ adduct)','Negative ionization (correct for H+ loss)')),
        h3("Max error (Da, 0 Da for theoritical masses):"),
        numericInput('tol', label='',0,min=0,max=0.5,step=0.001),
        h3("Enter Bibiserv2-Decomp job ID (Optional):"),
        textInput("DecompID", label='',value="Job ID")
),
    
    column(5,
        br(),
        h3("Faster peptide annotation"),
        em("For more than 500 masses"),
        br(),
        checkboxInput("checkbox", label = "Use DECOMP server (faster)", value = TRUE),
        br(),
        
        br(),
        tags$head(
          tags$style(HTML('#goButton{background-color:lightgreen}'))
        ),
        actionButton("goButton", "Start the run",style='padding:6px; font-size:150%'),
        br(),
        
        br(),
        tags$head(
          tags$style(HTML('#killButton{background-color:orange}'))
        ),
        actionButton("killButton", "Restart everything",style='padding:6px; font-size:150%'),
        br(),
        br(),
        br(),
        
        em('Messages from the server:'),
        br(),
        br(),
        textOutput("blank_message1")
    )),

    tabPanel("B) Mass decomposition",
        
        br(),
        
        tabsetPanel("",
        tabPanel("UAAC-Annotation",
          h3("Mass signals annotated to Unique Amino Acid Combination"),
          br(),
          dataTableOutput("table1")),
        
        tabPanel("MAAP-Annotation",
          h3("Mass signals annotated to Multiple Amino Acid Combination"),
          br(),
          dataTableOutput("table2"))),
        
        br(),
        downloadButton("downloadAnnotation", "Download decomposition results",style='padding:6px; font-size:150%')),
        
    tabPanel("C) Network construction",
        
        br(),
        column(8,
        h3("Please choose one or more options:"),
        
        HTML("<div style='height: 60px;'>"),
        checkboxGroupInput("visual", label = h6(""), choices = list("Include isolated nodes" = 1, "Remove triangles" = 2,"Remove free amino acids" = 3), selected = c(1,2)),
        HTML("</div>"),
        div(style="display: inline-block;vertical-align:middle; width: 300px;",checkboxGroupInput("visual2", label = h6(""), choices =list("Correlation between nodes" = 4), selected = 0)),
        div(style="display: inline-block;vertical-align:middle; width: 300px;",sliderInput("cor_min_max", label = h6(""), min = -1, max = 1, value = c(0.8,1),step=0.1)),
        br(),
        br(),
        htmlOutput("network_summary"),
        br(),
        dataTableOutput("table3"),
        br(),
        downloadButton("downloadNetwork", "Download network (edges) file",style='padding:6px; font-size:150%')),
        
        column(4,
        br(),
        br(),
        tags$head(
          tags$style(HTML('#goButtonbis{background-color:lightgreen}'))
        ),
        actionButton("goButtonbis", "Start network construction",style='padding:6px; font-size:150%'),
        br(),
        br(),
        tags$head(
          tags$style(HTML('#clearButton{background-color:orange}'))
        ),
        actionButton("clearButton", "Clear network",style='padding:6px; font-size:150%'),
        br(),
        br(),
        em('Messages from the server:'),
        textOutput("blank_message2"))),

    tabPanel("Network results", 
             
             # tags$style(type="text/css",
             #            ".shiny-output-error { visibility: hidden; }",
             #           ".shiny-output-error:before { visibility: hidden; }"
             # ),
             
        tabsetPanel("",
         tabPanel("Network visualization",
          column(3,
            h3("Peptide Degradation Network"),
            uiOutput("Shows1"),
            uiOutput("Shows2"),
            uiOutput("Shows3"),
            uiOutput("Coloring"),
            uiOutput("Color_type"),
            selectInput('Edge_Etiquette', 'Feature displayed when mouse over edges:',c("Amino acid(s) loss","Statistical correlation","Edge IDs")),
            br(),
            tags$head(
              tags$style(HTML('#goButton3{background-color:lightgreen}'))
            ),
            actionButton("goButton3", "Visualize entire network",style='padding:6px; font-size:150%')),
          
          column(9, 
            br(),
            visNetworkOutput("networkPlot",height = "1000px"))),
        
         tabPanel("Network analysis",
           br(),
           tags$head(
             tags$style(HTML('#goButton4{background-color:lightgreen}'))
           ),
           actionButton("goButton4", "Analyze the network",style='padding:6px; font-size:150%'),
           br(),
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
          br(),
          h3("Note: Same options will be used for visualizing subgraphs as for entire network !"),
          
          br(),
          tags$head(
              tags$style(HTML('#goButton5{background-color:lightgreen}'))
          ),
          actionButton("goButton5", "Extract high degree vertices",style='padding:6px; font-size:150%'),
          br(),
          br(),
          uiOutput("Controls"),
          selectInput('Order', 'Neighbourhood order:',1:5),
          visNetworkOutput("networkPlot3",height = "500px"),
          downloadButton("downloadNetwork3", "Download network (edges) file",style='padding:6px; font-size:150%'),
          br(),
          br(),
          br(),
          tags$head(
            tags$style(HTML('#goButton6{background-color:lightgreen}'))
          ),
          actionButton("goButton6", "Extract degradation chains",style='padding:6px; font-size:150%'),
          br(),
          br(),
          uiOutput("Chains"),
          visNetworkOutput("networkPlot1",height = "500px"),
          downloadButton("downloadNetwork1", "Download network (edges) file",style='padding:6px; font-size:150%')
      ))),
        
      tabPanel("Peptides in KEGG",
        br(),
        h3("Peptides in metabolic pathways"),
        selectInput('Organism', 'Choose your organism:',kegg_organism[,2]),
        uiOutput("Kegg_pathway"),
        htmlOutput('image')),
      
      tabPanel("Help",
              tabsetPanel("",
                    tabPanel("User Manual",
                        includeMarkdown("Help.Rmd")),
                    tabPanel("Troubleshooting",
                        includeMarkdown("Troubleshooting.Rmd")),
                    tabPanel("Example datasets",includeMarkdown("Examples.Rmd"))
      )),
      tabPanel("About",includeMarkdown("About.Rmd"))
      )
      
  
)


