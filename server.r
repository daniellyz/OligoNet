library(shiny)
options(warn=-1)
options(shiny.maxRequestSize=30*1024^2)
library(igraph)
library(KEGGREST)
require(visNetwork, quietly = TRUE)
require(stringr,  quietly = TRUE)
require(RCurl, quietly = TRUE)
require(DT, quietly = TRUE) 
#bibiserv2_2016-11-28_145427_MrRdE
#bibiserv2_2016-12-03_191118_7X2uK  # Without sugar/glutamine/asparagine
#bibiserv2_2016-12-05_161538_XtQxM # LC data 

source('script1.r')
source('script2.r')
source('script3.r')

kegg_organism=read.table("https://raw.githubusercontent.com/daniellyz/OligoNet/master/kegg_organisms.txt",sep="\t",header=F,stringsAsFactors = F)
kegg_cpds=read.table("https://raw.githubusercontent.com/daniellyz/OligoNet/master/kegg_cpds.txt",sep="\t",header=F,stringsAsFactors = F)
kegg_paths=read.table("https://raw.githubusercontent.com/daniellyz/OligoNet/master/kegg_pathway.txt",sep="\t",header=F,stringsAsFactors = F)

#deployApp(server="shinyapps.io",appName="OligoNet-V7")

#### Shiny ####
shinyServer(function(input, output) {
  
  selected_ex<-reactiveValues(a=0) # Which ex selected for execution, 0= no example
  annotated_data<-reactiveValues(k=list()) # If chosen from example, the decomposed data
  annotated_network<-reactiveValues(n=list()) # If chosen from example, the network
  
  monomers<-reactive({
  
    inFile3=input$file3
     
    if (is.null(inFile3)){
      monomers=read.csv('https://raw.githubusercontent.com/daniellyz/OligoNet/master/amino-acid-basic.txt',sep='\t',dec='.',header=F,stringsAsFactors = F)}
    else{monomers=read.csv(inFile3$datapath,sep='',dec='.',header=F,stringsAsFactors = F)}
    
    condensation=monomers[1,2]
    
    monomers=monomers[2:nrow(monomers),]

    monomers[,2]=monomers[,2]-condensation
    
    monomers=monomers[order(-monomers[,2]),] # Decreasing order
    
    elements=monomers[,1]
    
    list(tab=monomers,elements=elements,condensation=condensation)
  })
  
  find_annotation <- eventReactive(input$goButton,{

    annotated=NULL
    raw_data=NULL

    output_message=c()# Output messages for users
    
    inFile1=input$file1
    inFile2=input$file2
  
    if (!is.null(inFile1) && is.null(inFile2)) {inFile2=inFile1}
     
    if (!is.null(inFile1) && !is.null(inFile2)) { # Check the files users have uploaded
         
      raw_data=data.matrix(read.table(inFile1$datapath,sep='\t',dec=".",header=T,stringsAsFactors = F))
      additional_data=data.matrix(read.table(inFile2$datapath,sep='\t',dec=".",header=T,stringsAsFactors = F))
  
      if (!identical(raw_data[,1],additional_data[,1])){
        output_message=c(output_message,"The two data matrices must have identical ids")
        raw_data=addional_data=annotated=NULL}
      
      if (nrow(monomers()$tab)<2){
        output_message=c(output_message,"Not enough subunits for decomposition")
        raw_data=addional_data=annotated=NULL}
      
      if ((ncol(raw_data)<12) && (ncol(raw_data)>=3)){
        output_message=c(output_message,"Warning: At least 10 samples are needed to calculate correlations between nodes")}
      
      if (ncol(raw_data)<3){
        output_message=c(output_message,"File format not valid")
        raw_data=addional_data=annotated=NULL}
      }
    
    if (!is.null(raw_data) && !is.null(additional_data)){ # If everything is fine we proceed 
      
      if (ncol(raw_data)==3){raw_data=cbind(raw_data,III=rep(20,nrow(raw_data)))}
      
      dplace_raw_data=decimalnumcount(as.character(raw_data[2,2])) # decimal place
      dplace_monomers=decimalnumcount(as.character(monomers()$tab[2,2]))
      dplace=min(dplace_raw_data,dplace_monomers) # decimal place taken for rounding

      mass_list=round(raw_data[,2]-monomers()$condensation,dplace)
      
      tol1=0.01 # Default tolerance for Decomp (Theoritical)
      tol2=0.00001 # Default tolerance for filtering (Theoritical-Computational error)

      if (input$TE==2){tol2=input$tol}
      
      if (input$checkbox){# Use DECOMP server
      
        results=send_curl(mass_list,monomers()$tab,tol1,input$DecompID,dplace) # bielefeld results
      
        if (results$p3=="No response from the server"){ # If no output from the server
            annotated=NULL
            output_message=c(output_message,"The decomposition failed, please submit your job manually at http://bibiserv.techfak.uni-bielefeld.de/decomp")}
        
        else { # If everything OK
            output_message=c(output_message, "Decomposition succeeded! please download decomposition results!")
            output_massage=c(output_message,paste0("You Job id on DECOMP server is: ",results$id))
            annotated=peptide_annotation(raw_data,additional_data,results$p3,tol2)}}
      
      else { # Not use DECOMP server 
        annotated=peptide_annotation_slow(mass_list,dplace,monomers()$tab,raw_data,additional_data,tol2)       
        output_massage="Decomposition succeeded! please download decomposition results!"}}
 list(annotated=annotated,output_message=output_message,raw_data=raw_data,dplace=dplace,tol2=tol2)
  })
    
  output$summary <- renderPrint({
    if(selected_ex$a==0){
    output=find_annotation()$output_message
    }
    else {
    output=annotated_data$k
    output=output$output_message}
    cat(paste0(output,collapse="\n"))
    })
  
  output$table1 <- renderDataTable({
    
    if(selected_ex$a==0){found=find_annotation()}
    else{found=annotated_data$k}
    
   # save(found,file="FT-annotation.Rdata")
    tol2=found$tol2
    found=found$annotated

    UAAC=cbind(found$unique[,1:2],Peptide=found$unique[,"Peptide"])
    masslist=round(as.numeric(found$unique[,2]),4)

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
     else{cpd_code="0"}
    annotated_cpd=c(annotated_cpd,cpd_code)}
    UAAC <- cbind(UAAC,createLink(annotated_cpd))
    colnames(UAAC)[4]="KEGG"
    UAAC=datatable(UAAC,escape=c(TRUE, TRUE,TRUE, FALSE))

    return(UAAC)
  })
  
  output$table2 <- renderDataTable({
  
    if(selected_ex$a==0){found=find_annotation()}
    else{found=annotated_data$k}
    
    tol2=found$tol2
    found=found$annotated
    
    valid=found$all[,"NBP"]>1
    
    doubled=found$all[valid,]
    MAAP=cbind(doubled[,1:2],Peptide=doubled[,"Peptide"],NBP=doubled[,"NBP"])
    
    masslist=round(as.numeric(doubled[,2]),4)

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
      else{cpd_code="0"}
      annotated_cpd=c(annotated_cpd,cpd_code)}
    MAAP <- cbind(MAAP,createLink(annotated_cpd))
    colnames(MAAP)[5]="KEGG"
    MAAP=datatable(MAAP,escape=c(TRUE, TRUE,TRUE, TRUE, FALSE))
    return(MAAP)})
  
  load_network <- eventReactive(input$goButtonbis,{
    
      if(selected_ex$a==0){found=find_annotation()}
      else{found=annotated_data$k}
    
      found=found$annotated
      
      cor_min=-1
      delete_triangle=0
      
      if ("4" %in% input$visual){cor_min=input$cor_min}
      if ("2" %in% input$visual){delete_triangle=1}
      
      index=which(colnames(found$unique_add)==input$Etiquette)
      network=network_generator(found$unique,found$unique_add,cor_min,index,monomers()$elements,delete_triangle)
      
      if ("3" %in% input$visual){ # remove amino acids in the example files
        network=filter_amino_acid(network,monomers()$elements)}
      network
      
      })
  
  output$networkPlot <- renderVisNetwork({

    choose_to_show="1" %in% input$visual
    if(selected_ex$a==0){whole_network=load_network()}
    else{whole_network=annotated_network$n}
    
   #save(whole_network,file="FT-network.Rdata")
    if (!is.null(whole_network) && choose_to_show){
      visNetwork(whole_network$nodes, whole_network$links)}
  })
  
  output$Distribution_degree<-renderPlot({
    if(selected_ex$a==0){whole_network=load_network()}
    else{whole_network=annotated_network$n}
    
    deg_in <- degree(whole_network$g, mode="all")
    deg_in=deg_in[which(deg_in>=1)]
    plot(rank(-deg_in),as.numeric(deg_in),xlab="Rank",ylab="Degree")
})
  
  output$Distribution_edge<-renderPlot({
    if(selected_ex$a==0){whole_network=load_network()}
    else{whole_network=annotated_network$n}
    
    count_diff=table(whole_network$diff_list)
    plot(rank(-count_diff),as.numeric(count_diff),xlab="Rank",ylab="Occurences")
})  
  
  output$Top_edge<-renderPlot({
    if(selected_ex$a==0){whole_network=load_network()}
    else{whole_network=annotated_network$n}
    
    count_diff=table(whole_network$diff_list)
    count_diff=sort(count_diff,decreasing = T)
    most_rank=ceiling(length(count_diff)*0.2)
    barplot(count_diff[1:most_rank],ylab="Number of edges")
  })  
  
  output$Distribution_correlation<-renderPlot({
    if(selected_ex$a==0){whole_network=load_network()}
    else{whole_network=annotated_network$n}
    hist(whole_network$cor,xlab="Spearman correlation coefficient",ylab="Frequency",main="")
})  
  
  output$Distribution_length<-renderPlot({
    if(selected_ex$a==0){whole_network=load_network()}
    else{whole_network=annotated_network$n}
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
    barplot(table(dl-1),xlab="Path length between nodes",ylab="Frequency")
  })  
  
  find_cores<-reactive({
    
    if(selected_ex$a==0){whole_network=load_network()}
    else{whole_network=annotated_network$n}
    deg_in <- degree(whole_network$g, mode="in")
    #print(deg_in)
    w_in=which(deg_in>4)
    rw=match(w_in,whole_network$nodes$id)
    labels=whole_network$nodes$label[rw]
    list(labels=labels,w_in=w_in)
  })
  
  find_chains<-reactive({
    if(selected_ex$a==0){whole_network=load_network()}
    else{whole_network=annotated_network$n}
    all_chains=degrade_chain_all(whole_network,3)
    all_chains
  })
    
  output$networkPlot1 <- renderVisNetwork({

    if(selected_ex$a==0){whole_network=load_network()}
    else{whole_network=annotated_network$n}
    if (!is.null(whole_network)){
      chains_list=find_chains()
      paths_list=chains_list$paths
      wp=which(chains_list$start_labels==input$Start)
      if (length(wp)>0){
        path=paths_list[[wp]]
        chains=generate_degrade_chain(whole_network,path)
        visNetwork(chains$nodes, chains$links)}}
 })
  
 output$Controls <- renderUI({
    selectInput('Cores', 'Choose cores:',find_cores()$labels)
  })
 
 output$Chains <- renderUI({
   
   found_chains=find_chains()
   selectInput('Start', 'Choose chain that starts with:',found_chains$start_labels)
 })

 output$networkPlot3 <- renderVisNetwork({

   if(selected_ex$a==0){whole_network=load_network()}
   else{whole_network=annotated_network$n}
   if (!is.null(whole_network)){
    wp=which(find_cores()$labels==input$Cores)
    if (length(wp)>0){
      node_id=find_cores()$w_in[wp]
      cps=common_pattern_single(whole_network,node_id,'in',input$Order)
      visNetwork(cps$nodes, cps$links)}}
 })
 
 output$downloadAnnotation<-downloadHandler(
   
   filename = function() {
     "annotation.txt"},
   content = function(file) {
     if(selected_ex$a==0){found=find_annotation()}
     else{found=annotated_data$k}
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
 
 output$Shows <- renderUI({
   if(selected_ex$a==0){found=find_annotation()}
   else{found=annotated_data$k}
   found=found$annotated
   selectInput('Etiquette', 'Choose what to display when mouse on the nodes:',colnames(found$unique_add))
 })
 
 output$downloadNetwork<-downloadHandler(
   
   filename = function() {
     "edges.txt"},
   content = function(file) {
     if(selected_ex$a==0){whole_network=load_network()}
     else{whole_network=annotated_network$n}
     write.table(whole_network$edges, file,sep="\t",dec=",",row.names=F,col.names=T)}
 )
 
 output$downloadNetwork3<-downloadHandler(
   
   filename = function() {
     "high-degree.txt"},
   content = function(file) {
     if(selected_ex$a==0){whole_network=load_network()}
     else{whole_network=annotated_network$n}
     if (!is.null(whole_network)){
       wp=which(find_cores()$labels==input$Cores)
       node_id=find_cores()$w_in[wp]
       cps=common_pattern_single(whole_network,node_id,'all',input$Order)
       write.table(cps$edges,file,sep="\t",dec=",",row.names=F,col.names=T)}}
 )
 
 output$downloadNetwork1<-downloadHandler(
   
   filename = function() {
     "chain.txt"},
   content = function(file) {
     if(selected_ex$a==0){whole_network=load_network()}
     else{whole_network=annotated_network$n}
     if (!is.null(whole_network)){
       chains_list=find_chains()
       paths_list=chains_list$paths
       wp=which(chains_list$start_labels==input$Start)
       path=paths_list[[wp]]
       chains=generate_degrade_chain(whole_network,path)
       write.table(chains$edges,file,sep="\t",dec=",",row.names=F,col.names=T)}}
 )

output$Kegg_pathway <- renderUI({
  selectInput('pathway', 'Choose your pathway:',kegg_paths[,2])
})

output$image<-renderUI({
  
  ind=which(kegg_paths[,2]==input$pathway)
  code_path=kegg_paths[ind,1]
  code_path=str_extract(code_path,"[[:digit:]]+")
  wo=which(kegg_organism[,2]==input$Organism)
  organism=kegg_organism[wo,1] # Abbreviation of organism chosen
  code_path=paste0(organism,code_path)
  
  url=NULL
  if(selected_ex$a==0){found=find_annotation()}
  else{found=annotated_data$k}
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
  

observeEvent(input$Example1,{ 
    selected_ex$a=1
    load("FT-annotation.Rdata")
    load("FT-network.Rdata")
    annotated_data$k=found
    annotated_network$n=whole_network
})

observeEvent(input$Example2,{ 
  selected_ex$a=2
  load("LC-annotation.Rdata")
  load("LC-network.Rdata")
  annotated_data$k=found
  annotated_network$n=whole_network
})

observeEvent(input$clearButton,{ 
    selected_ex$a=0
    annotated_data$k=list()
    annotated_network$n=list()
})

})






