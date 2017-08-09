library(shiny)
library("V8")
library(shinyjs)
options(warn=-1)
options(shiny.maxRequestSize=60*1024^2)
library(igraph)
library(Hmisc)
library(KEGGREST)
require(visNetwork, quietly = TRUE)
require(stringr,  quietly = TRUE)
require(RCurl, quietly = TRUE)
require(DT, quietly = TRUE) 
#session$allowReconnect(TRUE)

#bibiserv2_2017-08-02_104105_jXVW5 # FT
#bibiserv2_2017-08-05_173101_ML09G # Wine LC
#bibiserv2_2017-08-02_152523_pr1kF # LC

source('script1.r')
source('script2.r')
source('script3.r')

kegg_organism=read.table("https://raw.githubusercontent.com/daniellyz/OligoNet/master/kegg_organisms.txt",sep="\t",header=F,stringsAsFactors = F)
kegg_cpds=read.table("https://raw.githubusercontent.com/daniellyz/OligoNet/master/kegg_cpds.txt",sep="\t",header=F,stringsAsFactors = F)
kegg_paths=read.table("https://raw.githubusercontent.com/daniellyz/OligoNet/master/kegg_pathway.txt",sep="\t",header=F,stringsAsFactors = F)

#deployApp(server="shinyapps.io",appName="OligoNet")

#### Shiny ####
shinyServer(function(input, output,clientData, session) {
  
monomers<-reactive({
  
    inFile3=input$file3
     
    if (is.null(inFile3)){
      monomers=read.csv('https://raw.githubusercontent.com/daniellyz/OligoNet/master/amino-acid-complete.txt',sep='\t',dec='.',header=F,stringsAsFactors = F)}
    else{monomers=read.csv(inFile3$datapath,sep='',dec='.',header=F,stringsAsFactors = F)}
    
    condensation=monomers[1,2]
    
    monomers=monomers[2:nrow(monomers),]

    monomers[,2]=monomers[,2]-condensation
    
    monomers=monomers[order(-monomers[,2]),] # Decreasing order
    
    elements=monomers[,1]
    
    list(tab=monomers,elements=elements,condensation=condensation)
  })
  
  
check_annotation <-eventReactive(input$goButton,{
  
    annotated=NULL
    raw_data=NULL
    additional_data=NULL
    dplace=NULL
    tol2=NULL
    mms=NULL
    valid=0 # If the data format is valid
      
    inFile1=input$file1
    inFile2=input$file2
    
    if ((is.null(inFile1)) && (input$blank_file1=="")){mms="Please load File 1!"}
    
    if (!is.null(inFile1)){ # Check file 1
      raw_data=read.csv(inFile1$datapath,sep='\t',dec=".",header=T,stringsAsFactors = F)}
    
    if ((is.null(inFile1)) && (input$blank_file1!="")){ # Read from the blank area instead
      input_str=input$blank_file1
      input_str=strsplit(input_str,"\n")
      input_str=input_str[[1]]
      header=strsplit(input_str[1],"\t")[[1]]
      ncolumns=length(header) # Number of columns of the matrix
      raw_data=sapply(input_str[2:length(input_str)],function(x) strsplit(x,"\t"))
      raw_data=do.call(rbind,raw_data)
      rownames(raw_data)=NULL
      if (nrow(raw_data)>1){
        I_data=apply(raw_data[,2:ncol(raw_data)],2,as.numeric)
        raw_data=data.frame(raw_data[,1],I_data,stringsAsFactors = F)}
      else{
        I_data=matrix(as.numeric(raw_data[,2:ncol(raw_data)]),nrow=1)
        raw_data=data.frame(raw_data[,1],I_data,stringsAsFactors = F)}
      colnames(raw_data)=header
    }
    
    if (!is.null(inFile2)){ # Check file 2
      additional_data=read.csv(inFile2$datapath,sep='\t',dec=".",header=T,stringsAsFactors = F)}
    
    if (!is.null(raw_data) && is.null(additional_data)){
      additional_data=raw_data}  # Default value for second file
    
    if (!is.null(raw_data) && !is.null(additional_data)){ # Start parameter checking
     
      if (ncol(raw_data)<3 || nrow(raw_data)<1){
        mms="File 1 must contain at least three columns and one row!"
        raw_data=addional_data=annotated=NULL}
      
      if (length(unique(raw_data[,1]))<length(raw_data[,1])){
        mms="Your data file must contain unique IDs!"
        raw_data=addional_data=annotated=NULL}
      
      if (!all.is.numeric(raw_data[,2])){
        mms="The second column representing mass values must contain numeric values!"
        raw_data=addional_data=annotated=NULL}
      
      if (!(all(raw_data[,2]>0 & raw_data[,2]<5000))){
        mms="Mass should be between 0 and 5000"
        raw_data=addional_data=annotated=NULL}
  
      if (nrow(monomers()$tab)<2){
        mms="Not enough subunits for decomposition"
        raw_data=addional_data=annotated=NULL}
        
      if (!(input$tol>=0 & input$tol<1)){
        mms="Mass error must between 0 and 1 Da"
        raw_data=addional_data=annotated=NULL}
    } 
    if (!is.null(raw_data) && !is.null(additional_data)){
      if (!identical(raw_data[,1],additional_data[,1])){
        mms="File 1 and File 2 must have identical ids"
        raw_data=addional_data=annotated=NULL}
      else{
      mms="File format valid!"
      valid=1
      }}
    list(raw_data=raw_data,additional_data=additional_data,message=mms,valid=valid)
  })
  
  
find_annotation <- eventReactive(input$goButton,{

    raw_data=check_annotation()$raw_data
    additional_data=check_annotation()$additional_data
    annotated=NULL
    dplace=NULL
    tol2=NULL
    mms=""

    if (!is.null(raw_data) && !is.null(additional_data)){ # If everything is fine we proceed 
      
      if ((ncol(raw_data)<6) && (ncol(raw_data)>=3)){
        mms="Warning: at least 4 samples are needed to calculate a meaningful statistical correlation"}
  
      if (as.numeric(raw_data[1,2])==round(as.numeric(raw_data[1,2]))){dplace=0} # integer
      else{
      dplace_raw_data=decimalnumcount(as.character(raw_data[1,2])) # decimal place
      dplace_monomers=decimalnumcount(as.character(monomers()$tab[2,2]))
      dplace=min(dplace_raw_data,dplace_monomers) # decimal place taken for rounding
      }  
      
      if (input$Scan=='Positive ionization (correct for H+ adduct)') {raw_data[,2]=raw_data[,2]-1.007276}
      if (input$Scan=='Negative ionization (correct for H+ loss)') {raw_data[,2]=raw_data[,2]+1.007276}
      
      mass_list=round(raw_data[,2]-monomers()$condensation,dplace)

      tol=input$tol # Tolerance 
      
        if (tol==0){
         if (dplace==0){tol=1}
        else{tol=10^(-dplace+1)}} # Default theoritical-computational error
      # 
      # if (dplace==0){tol=tol+1}
      # else {tol=tol+10^(-dplace+1)}
      
      monomers=monomers()$tab
      monomers[,2]=round(monomers[,2],dplace)
      condensation=monomers()$condensation
      condensation=round(condensation,dplace)
      
      if (input$checkbox){# Use DECOMP server
      
        results=send_curl(mass_list,monomers,tol,input$DecompID) # bielefeld results
      
        if (results$p3=="No response from the server"){ # If no output from the server
            annotated=NULL
            mms="The decomposition failed, but you can try to submit the job manually at http://bibiserv.techfak.uni-bielefeld.de/decomp.
            If your data has already been studied, you could try with its corresponding job id. Otherwise, you could
            run the job without using the DECOMP server (only recommended if you have fewer than 500 masses)"
          }
        
        else { # If everything OK
          mms="Decomposition succeeded! Please check the results on the next tab-panel! Your job id is: "
          mms=paste0(mms,results$id,", you could enter this id if you want to re-run your data.")
          annotated=peptide_annotation(raw_data,additional_data,results$p3,tol,monomers,condensation)
          }}
      
      else { # Not use DECOMP server 
        annotated=peptide_annotation_slow(mass_list,monomers,raw_data,additional_data,tol,condensation)       
        mms="Decomposition succeeded! Please check the results on the next tab-panel!"
        }
  }
 list(annotated=annotated,raw_data=raw_data,dplace=dplace,tol2=tol,message=mms)
  })
  
observeEvent(input$goButton,{

  withProgress({
    setProgress(message="Check data format...")
    Sys.sleep(1)
    setProgress(message=check_annotation()$message)
    Sys.sleep(1)
    if (check_annotation()$v==1){
      setProgress(message="Decomposition started...")
      Sys.sleep(1)
      setProgress(message=find_annotation()$message)
    }
  })
  if (check_annotation()$v==1){
  output$blank_message1<-renderText({find_annotation()$message})}      
  if (check_annotation()$v==0){
    updateActionButton(session, "goButton",label = "Try again!")
    output$blank_message1<-renderText({check_annotation()$message})
  }
})

observeEvent(input$killButton,{shinyjs::js$refresh()})

output$table1 <- renderDataTable({
    
  withProgress({
    setProgress(message="Generating annotation results...")
    Sys.sleep(1)
    found=find_annotation()
    
    tol2=found$tol2
    found=found$annotated
    
    if (nrow(found$unique)==0){UAAC=NULL}
    
    if (nrow(found$unique)>1){
     UAAC=cbind(found$unique[,1:2],Peptide=found$unique[,"Peptide"],Error_ppm=found$unique[,"PPM"])}

    if (nrow(found$unique)==1){
       UAAC=c(found$unique[,1:2],found$unique[,"Peptide"],found$unique[,"PPM"])   
       UAAC=matrix(UAAC,nrow=1,dimnames=list(NULL,c("ID","Mass","Peptide","Error_ppm")))}
    
    UAAC=datatable(UAAC,escape=c(TRUE,TRUE,TRUE,TRUE), rownames = F)

    return(UAAC)
  })
  })
  
output$table2 <- renderDataTable({
  
  withProgress({
    setProgress(message="Generating annotation results...")
    Sys.sleep(1)
  
    found=find_annotation()
    
    tol2=found$tol2
    found=found$annotated
    
    valid=found$all[,"NBP"]>1
    
    doubled=found$all[valid,]
    
    if (nrow(doubled)==0){MAAP=NULL}
  
    if (nrow(doubled)>1){
    MAAP=cbind(doubled[,1:2],Peptide=doubled[,"Peptide"],NBP=doubled[,"NBP"],Error_ppm=doubled[,"PPM"])}
    
    if (nrow(doubled)==1){
      MAAP=c(doubled[,1:2],doubled[,"Peptide"],doubled[,"NBP"],doubled[,"PPM"])   
      MAAP=matrix(MAAP,nrow=1,dimnames=list(NULL,c("ID","Mass","Peptide","NBP","Error_ppm")))}
    
    MAAP=datatable(MAAP,escape=c(TRUE,TRUE,TRUE,TRUE,TRUE), rownames = F)
    return(MAAP)})
})
  
observeEvent(input$goButtonbis ,{

  withProgress({
  found=find_annotation()
  found=found$annotated
  found_unique=found$unique
  
  if (!is.null(found_unique)){ 
    if (nrow(found$unique)<2){
        setProgress(message="At least 2 UAAC signals required for network construction...")}
    else{
      setProgress(message="Start network construction...")
      whole_network=load_network()
      save(whole_network,file="FT-modified.rdata")
      Sys.sleep(1)
      if (length(whole_network$cor)==0){
        setProgress(message="No network structure found in your data!")}
      else{
        setProgress(message="Network construction succeeded!")
      }}}
  Sys.sleep(1)
  setProgress(message="Generating edge table...")
  Sys.sleep(1)
 })
 if (is.null(found_unique)){
    output$blank_message2<-renderText("Please start the run first in tab-panel A)")}  
 else{
  if (nrow(found$unique)<2){
    output$blank_message2<-renderText("At least 2 UAAC signals required for network construction...")}      
  else{
    if (length(whole_network$cor)==0){
      output$blank_message2<-renderText("No network structure found in your data!")}
    else{
     output$blank_message2<-renderText("Network construction succeeded! You can check and download the network file below for visualization in Gephi/Cytoscape, or go to the next tab-panel to visualize or analyze the network!")
     
     N1=length(whole_network$nodes$id)
     N2=length(whole_network$cor)
     m1=paste("<b>","Summary of the network:","</b>")
     m2=paste0("Your network contains ",N1," nodes and ",N2," edges.")
     if (ncol(found$unique)<9){m3="Statistical correlation might not be meaningful since you have fewer than 4 samples!"}
     else{m3=""}
     m4="List of edges in your network:"
     output$network_summary=renderUI({HTML(paste(m1,m2,m3,m4,sep="<br/>"))})
     #random_edges=sample(1:nrow(whole_network$edges),10)
     
     EDGES=whole_network$edges
     EDGES=datatable(EDGES,escape=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE))
     output$table3 <- renderDataTable({EDGES})
    }}}
})

observeEvent(input$clearButton,{
  output$network_summary<-renderUI({HTML(NULL)})
  output$blank_message2<-renderText("")
  output$table3<-renderDataTable(NULL)
  output$networkPlot <- renderVisNetwork(NULL)
  output$Distribution_degree<-renderPlot(NULL)
  output$Distribution_edge<-renderPlot(NULL)
  output$Top_edge<-renderPlot(NULL)
  output$Distribution_correlation<-renderPlot(NULL)
  output$Distribution_length<-renderPlot(NULL)
 # output$networkPlot1 <- renderVisNetwork(NULL)
  output$Controls <- renderUI(NULL)
  output$Chains <- renderUI(NULL)
#  output$networkPlot3 <- renderVisNetwork(NULL)
})

load_network <- eventReactive(input$goButtonbis,{
  output$network_summary<-renderUI({HTML(NULL)})
  output$blank_message2<-renderText("")
  output$table3<-renderDataTable(NULL)
  output$networkPlot <- renderVisNetwork(NULL)
  output$Distribution_degree<-renderPlot(NULL)
  output$Distribution_edge<-renderPlot(NULL)
  output$Top_edge<-renderPlot(NULL)
  output$Distribution_correlation<-renderPlot(NULL)
  output$Distribution_length<-renderPlot(NULL)
#  output$networkPlot1 <- renderVisNetwork(NULL)
  output$Controls <- renderUI(NULL)
  output$Chains <- renderUI(NULL)
#  output$networkPlot3 <- renderVisNetwork(NULL)  
    
  found=find_annotation()
  found=found$annotated
  network=NULL
    
  # if (!all.is.numeric(found$unique[,1])){
  #   found$unique[,1]=1:nrow(found$unique)
  #   found$unique_add[,1]=1:nrow(found$unique) # Replace by numeric
  # }
  
  
  if (nrow(found$unique)>1){
    cor_min=-100
    cor_max=100
      
    delete_triangle=0
    isolated=0
    filter_aa=0
    
    if ("4" %in% input$visual2){
      cor_min=input$cor_min_max[1]
      cor_max=input$cor_min_max[2]}
    
    if ("2" %in% input$visual){delete_triangle=1}
    if ("1" %in% input$visual){isolated=1}
    if ("3" %in% input$visual){filter_aa=1}

    #network_generator<-function(unique,unique_add,cor_min,cor_max,elements,delete_triangle,isolated){
      
    network=network_generator(found$unique,found$unique_add,cor_min,cor_max,monomers()$elements,delete_triangle,isolated,filter_aa)
    }
    network
})


output$Shows1 <- renderUI({
  names_list=NULL
  found=find_annotation()
  found=found$annotated
  if (!is.null(found)){
    names_list=unique(c(colnames(found$unique),colnames(found$unique_add)))}
  names_list=c("Peptide",names_list)
  selectInput('Display', 'Feature displayed next to each node:',names_list)
})

output$Shows2 <- renderUI({
  names_list=NULL
  found=find_annotation()
  found=found$annotated
  if (!is.null(found)){
    names_list=unique(c(colnames(found$unique),colnames(found$unique_add)))}
  names_list=c("<None>",names_list)
  selectInput('Etiquette', 'Feature displayed when mouse over nodes:',names_list)
})

output$Shows3 <- renderUI({
  names_list=NULL
  found=find_annotation()
  found=found$annotated
  if (!is.null(found)){
    names_list=unique(c(colnames(found$unique),colnames(found$unique_add)))}
  names_list=c("<None>",names_list)
  selectInput('NetworkColor', 'Feature used for node coloring:',names_list)
})

output$Color_type <- renderUI({
  
  names_list=NULL
  found=find_annotation()
  whole_network=load_network()
  
  found=found$annotated
  found=cbind(found$unique,found$unique_add)
  found=found[whole_network$matched,]

  index2=which(colnames(found)==input$NetworkColor)[1]
  if (!is.na(index2)){
  variable_to_color=found[,index2]
  variable2=unique(variable_to_color)
  
  if (length(variable2)>1){
    if (!all.is.numeric(c(variable2))){
      names_list=c()
      for (i in 1:length(variable2)){
        names=paste0("Color the group ",variable2[i]," red")
        names_list=c(names_list,names)}}
    else{
      names_list="Default coloring: quantitaive values"
      if (length(variable2)<6){
        for (i in 1:length(variable2)){
          names=paste0("Color the group ",variable2[i]," red")
          names_list=c(names_list,names)}}  }}}
  selectInput('Modality', 'Type of node coloring:',names_list)
})

network_modifier <- reactive({

    whole_network=load_network()
    if (!is.null(whole_network)){
    
    found=find_annotation()
    
    found=found$annotated
    found=cbind(found$unique,found$unique_add)
    ranked=order(found[,2])
    found=found[ranked,]
    found=found[whole_network$matched,]
  
    index0=which(colnames(found)==input$Display)[1]
    index1=which(colnames(found)==input$Etiquette)[1]
    index2=which(colnames(found)==input$NetworkColor)[1]
    
    if (!is.na(index0)){whole_network$nodes$label <- found[,index0]}
    if (!is.na(index1)){whole_network$nodes$title <- found[,index1]}
    else{whole_network$nodes$title<-NA}
  
    col_list=rep("blue",nrow(found))
    if (!is.na(index2)){
    variable_to_color=found[,index2]
    variable2=unique(variable_to_color)
    
    if (length(variable2)>1){
      if (!all.is.numeric(c(variable2))){
        names_list=c()
        for (i in 1:length(variable2)){
          names=paste0("Color the group ",variable2[i]," red")
          names_list=c(names_list,names)}}
      else{
        names_list="Default coloring: quantitaive values"
        if (length(variable2)<6){
          for (i in 1:length(variable2)){
            names=paste0("Color the group ",variable2[i]," red")
            names_list=c(names_list,names)}}}
      
    w=which(names_list==input$Modality)
    if (names_list[w]=="Default coloring: quantitaive values"){
      variable_to_color=rank(as.numeric(variable_to_color),ties.method="first")
      colfunc <- colorRampPalette(c("lightblue", "darkblue"))
      color_levels=colfunc(max(variable_to_color)) 
      col_list=color_levels[variable_to_color]}
    else{
      #if (all.is.numeric(c(variable2))){w=w-1} # Quatitative values
      var=variable2[w-1]
      group_selected=which(variable_to_color==var)
      col_list[group_selected]="red"}}}
    
    whole_network$nodes$color.background=col_list
    
    ee=c("Amino acid(s) loss","Statistical correlation","Edge IDs")
    if (input$Edge_Etiquette==ee[2]){whole_network$links$title=whole_network$edges[,7]}
    if (input$Edge_Etiquette==ee[3]){whole_network$links$title=whole_network$edges[,1]}}
    
    whole_network
})

observeEvent(input$goButton3,{
  withProgress({
   setProgress(message="Generate the entire network...")
   Sys.sleep(1)
   whole_network=network_modifier()
   if (!is.null(whole_network)){
    output$networkPlot<- renderVisNetwork({visNetwork(whole_network$nodes, whole_network$links)})
  }
   setProgress(message="The network will soon show up in your browser...")
   Sys.sleep(1)
   setProgress(message="Please be patient! This can be long depending on your browser.")
   Sys.sleep(5)
  })
})

observeEvent(input$goButton4,{
  withProgress({
  setProgress(message="Analyzing global network topology...")
  Sys.sleep(2)
    
  whole_network=load_network()

  deg_in <- degree(whole_network$g, mode="all")
  deg_in=deg_in[which(deg_in>=1)]
  output$Distribution_degree<-renderPlot({plot(rank(-deg_in),as.numeric(deg_in),xlab="Rank",ylab="Degree")})
  
  count_diff=table(whole_network$diff_list)
  output$Distribution_edge<-renderPlot({plot(rank(-count_diff),as.numeric(count_diff),xlab="Rank",ylab="Occurences")})
  
  count_diff=sort(count_diff,decreasing = T)
  most_rank=ceiling(length(count_diff)*0.2)
  output$Top_edge<-renderPlot({barplot(count_diff[1:most_rank],ylab="Number of edges")})  
  
  setProgress(message="Analyzing statistical properties...")
  Sys.sleep(2)
  output$Distribution_correlation<-renderPlot({hist(as.numeric(whole_network$cor),xlab="Spearman correlation coefficient",ylab="Frequency",main="")})
  
  setProgress(message="Analyzing path lengths in the network...")
  g=whole_network$g
  nodes=whole_network$nodes
  dl=c()
  for (n in 1:(nrow(nodes)-1)){
    for (m in ((n+1):nrow(nodes))){
      checked_nodes=try(all_simple_paths(g,from=nodes[m,1],to=nodes[n,1]),silent=T)
      if ((length(checked_nodes)>0) & (class(checked_nodes)!="try-error")){
        dl=c(dl,sapply(checked_nodes,length))}
    }
  }
  output$Distribution_length<-renderPlot({barplot(table(dl-1),xlab="Path length between nodes",ylab="Frequency")}) 
  
  setProgress(message="Network analysis succeeded!")
  Sys.sleep(1)})
  
})

observeEvent(input$goButton5,{
  
  withProgress({
  whole_network0=load_network()
  if (!is.null(whole_network0)){
    setProgress(message="Looking for high degree nodes in the network...")
    Sys.sleep(1)
    deg_in <- degree(whole_network0$g, mode="in")
    w_in=which(deg_in>4)
    rw=match(names(w_in),whole_network0$nodes$id)
    labels=whole_network0$nodes$label[rw]
    output$Controls <- renderUI({selectInput('Cores', 'Choose center nodes:',labels)})}
  if (length(labels)==0){
    setProgress(message="No nodes in the network have a degree >4 !")
    Sys.sleep(1)}
  })
})

output$networkPlot3 <- renderVisNetwork({
  whole_network0=load_network()
  if (!is.null(whole_network0)){
    deg_in <- degree(whole_network0$g, mode="in")
    w_in=which(deg_in>4)
    rw=match(names(w_in),whole_network0$nodes$id)
    labels=whole_network0$nodes$label[rw]
    wp=which(labels==input$Cores)
    if (length(wp)>0){
      whole_network=network_modifier()
      save(whole_network,file="whole_network2.rdata")
      node_id=names(w_in)[wp]
      cps=common_pattern_single(whole_network,node_id,'in',input$Order)
      visNetwork(cps$nodes, cps$links)}}
})

output$downloadNetwork3<-downloadHandler(
  filename = function() {
    "high-degree.txt"},
  content = function(file) {
    whole_network0=load_network()
    if (!is.null(whole_network0)){
      deg_in <- degree(whole_network0$g, mode="in")
      w_in=which(deg_in>4)
      rw=match(names(w_in),whole_network0$nodes$id)
      labels=whole_network0$nodes$label[rw]
      wp=which(labels==input$Cores)
      if (length(wp)>0){
        whole_network=network_modifier()
        node_id=names(w_in)[wp]
        cps=common_pattern_single(whole_network,node_id,'in',input$Order)
        write.table(cps$edges,file,sep="\t",dec=",",row.names=F,col.names=T)}}}
)

find_chains <- eventReactive(input$goButton6,{
  whole_network0=load_network()
  all_chains=NULL
  withProgress({
  setProgress(message="Looking for long paths (length>3) in the network...")
  Sys.sleep(1)
  if (!is.null(whole_network0)){
  all_chains=degrade_chain_all(whole_network0,3)
  setProgress(message="All long paths are found! Please choose which chain you want to visualize!")}
  else{setProgress(message="No long paths detected!")}
  Sys.sleep(2)
  })
  all_chains
})
 
observeEvent(input$goButton6,{
  if (!is.null(find_chains())){
  output$Chains <- renderUI({selectInput('Start', 'Choose chain that you want to visualize:',find_chains()$start_labels)})
  }
})

observeEvent(input$Start,{
  whole_network0=load_network()
  if (!is.null(whole_network0)){
    paths_list=find_chains()$paths
    wp=which(find_chains()$start_labels==input$Start)
    if (length(wp)>0){
        whole_network=network_modifier()
        path=paths_list[[wp]]
        chains=generate_degrade_chain(whole_network,path)
        output$networkPlot1 <- renderVisNetwork({visNetwork(chains$nodes, chains$links)})
    }}
})

output$downloadNetwork1<-downloadHandler(
  
  filename = function() {
    "chain.txt"},
  content = function(file) {
    whole_network0=load_network()
    if (!is.null(whole_network0)){
      paths_list=find_chains()$paths
      wp=which(find_chains()$start_labels==input$Start)
      if (length(wp)>0){
        whole_network=network_modifier()
        path=paths_list[[wp]]
        chains=generate_degrade_chain(whole_network,path)
        write.table(chains$edges,file,sep="\t",dec=",",row.names=F,col.names=T)}}}
)

output$downloadAnnotation<-downloadHandler(
   filename = function() {
     "annotation.txt"},
   content = function(file) {
     found=find_annotation()
     tol2=found$tol2
     found=found$annotated
     masslist=round(as.numeric(found$all[,2]),digits=4) # only masses annotated to unique combination..
     annotated_index=list()
     for (m in 1:length(masslist)){
       valid=which(abs(kegg_cpds[,2]-masslist[m])<=tol2)
       annotated_index[[m]]=valid}
     annotated_cpd=c()
     for (i in 1:length(annotated_index)){
       if (length(annotated_index[[i]])==1){
         cpd_code=kegg_cpds[annotated_index[[i]],1]
         cpd_code=strsplit(cpd_code,"cpd:")[[1]][2]}
       else if (length(annotated_index[[i]])>1){
         cpd_code=kegg_cpds[annotated_index[[i]],1]
         cpd_code=sapply(cpd_code,function(x) strsplit(x,"cpd:")[[1]][2])
         cpd_code=paste0(cpd_code,collapse = " , ")}
       else{cpd_code="Non metabolites found"}
       annotated_cpd=c(annotated_cpd,cpd_code)}
  all_annotated= cbind(found$all_add,KEGG=annotated_cpd)     
  write.table(all_annotated, file,sep="\t",dec=",",row.names=F,col.names=T)}
)

output$downloadNetwork<-downloadHandler(
   
   filename = function() {
     "edges.txt"},
   content = function(file) {
     whole_network=load_network()
     write.table(whole_network$edges, file,sep="\t",dec=",",row.names=F,col.names=T)}
)
 
output$Kegg_pathway <- renderUI({
  selectInput('pathway', 'Choose your pathway*:',kegg_paths[,2])
})

output$image<-renderUI({
  
  ind=which(kegg_paths[,2]==input$pathway)
  code_path=kegg_paths[ind,1]
  code_path=str_extract(code_path,"[[:digit:]]+")
  wo=which(kegg_organism[,2]==input$Organism)
  organism=kegg_organism[wo,1] # Abbreviation of organism chosen
  code_path=paste0(organism,code_path)
  
  url=NULL
  found=find_annotation()
  tol2=found$tol2
  found=found$annotated
  raw_data=found$all
  raw_data_peptide=raw_data[which(raw_data[,"NBP"]>0),] # Peptide features
  ppmass=round(as.numeric(raw_data_peptide[,2]),4)
  raw_data_non_peptide=raw_data[which(raw_data[,"NBP"]==0),]
  npmass=round(as.numeric(raw_data_non_peptide[,2]),4)

  ll=try(keggGet(code_path),silent=T)
  if (class(ll)!="try-error"){
      cpds=try(names(ll[[1]]$COMPOUND),silent=T)
      if ((class(cpds)!="NULL") && (class(cpds)!="try-error")) {
        cpds=paste0('cpd:',cpds)
        matched_cpds_ind=match(cpds,kegg_cpds[,1])
        matched_cpds_ind=matched_cpds_ind[!is.na(matched_cpds_ind)]
        cpds_mass_list=round(as.numeric(kegg_cpds[matched_cpds_ind,2]),4)
        cpds_list=kegg_cpds[matched_cpds_ind,1]
        
        nodes_peptide=c()
        nodes_non_peptide=c()
        for (k in 1:length(cpds_list)){
          vp=which(abs(cpds_mass_list[k]-ppmass)<=tol2)
          if (length(vp>0)){nodes_peptide=c(nodes_peptide,cpds_list[k])}
          vn=which(abs(cpds_mass_list[k]-npmass)<=tol2)
          if (length(vn>0)){nodes_non_peptide=c(nodes_non_peptide,cpds_list[k])}}
        
        if ((length(nodes_peptide)>0) || length(nodes_non_peptide)>0){
        col1=rep('red',length(nodes_peptide))
        col2=rep('navy',length(nodes_non_peptide))
        url= color.pathway.by.objects(code_path,c(nodes_peptide,nodes_non_peptide),c(col1,col2), c(col1,col2))}
      }
    }
  tags$img(src=url)
  })
  
})






