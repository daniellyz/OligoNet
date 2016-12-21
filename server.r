library(shiny)
options(warn=-1)
options(shiny.maxRequestSize=30*1024^2)
library(igraph)
library(KEGGREST)
require(visNetwork, quietly = TRUE)
require(stringr,  quietly = TRUE)
require(RCurl, quietly = TRUE)
library(DT) 
#bibiserv2_2016-11-28_145427_MrRdE
#bibiserv2_2016-12-03_191118_7X2uK  # Without sugar/glutamine/asparagine
#bibiserv2_2016-12-05_161538_XtQxM # LC data 

source('script1.r')
source('script2.r')
kegg_organism=read.table("kegg_organisms.txt",sep="\t",header=F,stringsAsFactors = F)
kegg_cpds=read.table("kegg_cpds.txt",sep="\t",header=F,stringsAsFactors = F)

#deployApp(server="shinyapps.io",appName="OligoNet-V7")

#### Shiny ####
shinyServer(function(input, output) {
  
  
  monomers<-reactive({
  
    inFile3=input$file3
     
    if (is.null(inFile3)){
      monomers=read.csv('https://github.com/daniellyz/OligoNet/blob/master/amino-acid-basic.txt',sep='\t',dec=',',header=F,stringsAsFactors = F)}
    else{monomers=read.csv(inFile3$datapath,sep='',dec=',',header=F,stringsAsFactors = F)}
    
    condensation=monomers[1,2]
    
    monomers=monomers[2:nrow(monomers),]
    
    monomers[,2]=monomers[,2]-condensation
    
    elements=monomers[,1]
    
    list(tab=monomers,elements=elements,condensation=condensation)
  })
  
  find_annotation <- eventReactive(input$goButton,{

    annotated=NULL
    raw_data=NULL
    
    inFile1=input$file1
    inFile2=input$file2
    
    output_message=c()# Output messages for users
    
    if (!is.null(inFile1) && is.null(inFile2)) {inFile2=inFile1}
     
    if (!is.null(inFile1) && !is.null(inFile2)) { # Check the files users have uploaded
         
      raw_data=data.matrix(read.table(inFile1$datapath,sep='\t',dec=',',header=T,stringsAsFactors = F))
      additional_data=data.matrix(read.table(inFile2$datapath,sep='\t',dec=',',header=T,stringsAsFactors = F))
  
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
      
      dplace=decimalnumcount(as.character(raw_data[1,2])) # decimal place
        
      mass_list=raw_data[,2]-round(monomers()$condensation,dplace)
      
      tol1=0.01 # Default tolerance for Decomp (Theoritical)
      
      tol2=1e-5 # Default tolerance for filtering (Theoritical-Computational error)
      
      if (input$TE==2){tol2=input$tol}
      
      results=send_curl(mass_list,monomers()$tab,tol1,input$DecompID,dplace) # bielefeld results
      
      if (results=="No response from the server"){ # If no output from the server
          annotated=NULL
          output_message=c(output_message,"The decomposition failed, please submit your job manually at http://bibiserv.techfak.uni-bielefeld.de/decomp")}
        
      else { # If everything OK
          output_message=c(output_message, "Decomposition succeeded! please download decomposition results!")
          annotated=peptide_annotation(raw_data,additional_data,results,tol2)}
    }
 save(annotated,file="raw.Rdata")
 list(annotated=annotated,output_message=output_message,raw_data=raw_data,dplace=dplace,tol2=tol2)
  })
    
  
  output$summary <- renderPrint({
    output=find_annotation()$output_message
    cat(paste0(output,collapse="\n"))
    })
  
  output$table1 <- renderDataTable({

    found=find_annotation()$annotated
    dplace=find_annotation()$dplace
    
    UAAC=cbind(found$unique[,1:2],Peptide=found$unique[,"Peptide"])
    
    masslist=as.numeric(found$unique[,2])
    if (dplace>4){masslist=round(masslist,digits=4)} # only masses annotated to unique combination..
    if (dplace<4){kegg_cpds[,2]=round(kegg_cpds[,2],dplace)}
    
    annotated_index=lapply(masslist,function(x) which(kegg_cpds[,2]==x)) # index in kegg that fits measured mass
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
    
    found=find_annotation()$annotated
    dplace=find_annotation()$dplace
  
    valid=found$all[,"NBP"]>1
    doubled=found$all[valid,]
    MAAP=cbind(doubled[,1:2],Peptide=doubled[,"Peptide"],NBP=doubled[,"NBP"])
    
    masslist=as.numeric(doubled[,2])
    if (dplace>4){masslist=round(masslist,digits=4)} # only masses annotated to unique combination..
    if (dplace<4){kegg_cpds[,2]=round(kegg_cpds[,2],dplace)}
    
    annotated_index=lapply(masslist,function(x) which(kegg_cpds[,2]==x)) # index in kegg that fits measured mass
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
    
      cor_min=-1
      
      if ("3" %in% input$visual){cor_min=input$cor_min}
  
      found=find_annotation()$annotated
      
      index=which(colnames(found$unique_add)==input$Etiquette)
        
      network=network_generator(found$unique,found$unique_add,cor_min,index,monomers()$elements)
      
      if ("2" %in% input$visual){
        network=filter_amino_acid(network,monomers()$elements)}
      save(network,file="network_filtered.Rdata")
  network})
  
  
  output$Distribution_degree<-renderPlot({
    whole_network=load_network()
    deg_in <- degree(whole_network$g, mode="all")
    deg_in=deg_in[which(deg_in>=1)]
    plot(rank(-deg_in),as.numeric(deg_in),xlab="Rank",ylab="Degree")  
})
  
  output$Distribution_edge<-renderPlot({
    whole_network=load_network()
    count_diff=table(whole_network$diff_list)
    plot(rank(-count_diff),as.numeric(count_diff),xlab="Rank",ylab="Occurences")
})  
  
  output$Top_edge<-renderPlot({
    whole_network=load_network()
    count_diff=table(whole_network$diff_list)
    count_diff=sort(count_diff,decreasing = T)
    most_rank=ceiling(length(count_diff)*0.2)
    barplot(count_diff[1:most_rank],ylab="Number of edges")
  })  
  
  output$Distribution_correlation<-renderPlot({
    whole_network=load_network()
    hist(whole_network$cor,xlab="Spearman correlation coefficient",ylab="Frequency",main="")
})  
  
  output$Distribution_length<-renderPlot({
    whole_network=load_network()
    g=whole_network$g
    nodes=whole_network$nodes
    dl=c()
    for (n in 1:(nrow(nodes)-1)){
      for (m in ((n+1):nrow(nodes))){
        checked_nodes=all_simple_paths(g,from=nodes[m,1],to=nodes[n,1])
        if (length(checked_nodes)>0){
          dl=c(dl,sapply(checked_nodes,length))}
      }
    }
    barplot(table(dl-1),xlab="Path length between nodes",ylab="Frequency")
  })  
  
  find_cores<-reactive({
    
    whole_network=load_network()
    deg_in <- degree(whole_network$g, mode="all")
    #print(deg_in)
    w_in=which(deg_in>4)
    rw=match(w_in,whole_network$nodes$id)
    labels=whole_network$nodes$label[rw]
    list(labels=labels,w_in=w_in)
  })
  
  find_chains<-reactive({
    whole_network=load_network()
    all_chains=degrade_chain_all(whole_network,3)
    all_chains
  })
    
  
  output$networkPlot1 <- renderVisNetwork({

    whole_network=load_network()
    if (!is.null(whole_network)){
    chains_list=find_chains()
    paths_list=chains_list$paths
    wp=which(chains_list$start_labels==input$Start)
    path=paths_list[[wp]]
    chains=generate_degrade_chain(whole_network,path)
    visNetwork(chains$nodes, chains$links)}
 })
  
  
 output$Controls <- renderUI({
    selectInput('Cores', 'Choose cores:',find_cores()$labels)
  })
 
 output$Chains <- renderUI({
   
   found_chains=find_chains()
   selectInput('Start', 'Choose chain that starts with:',found_chains$start_labels)
 })

 output$networkPlot3 <- renderVisNetwork({

   whole_network=load_network()
   if (!is.null(whole_network)){
    wp=which(find_cores()$labels==input$Cores)
    node_id=find_cores()$w_in[wp]
    cps=common_pattern_single(whole_network,node_id,'all',input$Order)
    visNetwork(cps$nodes, cps$links)}
 })
 
 output$downloadAnnotation<-downloadHandler(
   
   filename = function() {
     "annotation.txt"},
   content = function(file) {
     found=find_annotation()$annotated
     masslist=round(as.numeric(found$all[,2]),digits=4) # only masses annotated to unique combination..
     annotated_index=lapply(masslist,function(x) which(kegg_cpds[,2]==x)) # index in kegg that fits measured mass
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
  all_annotated= cbind(found$all,KEGG=annotated_cpd)     
  write.table(all_annotated, file,sep="\t",dec=",",row.names=F,col.names=T)}
)
 
 output$Shows <- renderUI({
   found=find_annotation()$annotated
   selectInput('Etiquette', 'Choose what to display when mouse on the nodes:',colnames(found$unique_add))
 })
 
 output$downloadNetwork<-downloadHandler(
   
   filename = function() {
     "edges.txt"},
   content = function(file) {
     whole_network=load_network()
     write.table(whole_network$edges, file,sep="\t",dec=",",row.names=F,col.names=T)}
 )
 
 output$downloadNetwork3<-downloadHandler(
   
   filename = function() {
     "high-degree.txt"},
   content = function(file) {
     whole_network=load_network()
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
     whole_network=load_network()
     if (!is.null(whole_network)){
       chains_list=find_chains()
       paths_list=chains_list$paths
       wp=which(chains_list$start_labels==input$Start)
       path=paths_list[[wp]]
       chains=generate_degrade_chain(whole_network,path)
       write.table(chains$edges,file,sep="\t",dec=",",row.names=F,col.names=T)}}
 )

 output$networkPlot <- renderVisNetwork({
   choose_to_show="1" %in% input$visual
   whole_network=load_network()
    if (!is.null(whole_network) && choose_to_show){
      visNetwork(whole_network$nodes, whole_network$links)}
  })
 
load_KEGG<- eventReactive(input$goButton3,{  

   found=find_annotation()$annotated
   raw=find_annotation()$raw_data
   
   row_list_peptide=which(as.numeric(found$all[,"NBP"])>0)
   masslist=round(as.numeric(found$all[row_list_peptide,2]),digits=4) # only masses annotated to unique combination..

   index_cpd=match(masslist,kegg_cpds[,2])
   index_cpd=index_cpd[which(!is.na(index_cpd))]   
   code_cpd=kegg_cpds[index_cpd,1]
   
   pathways=sapply(code_cpd,function(x) keggGet(x)[[1]]$PATHWAY)
   pn=unique(unlist(lapply(pathways,function(x) names(x))))
   pathways=unique(unlist(pathways)) # Find all pathway involved
   
   code_path=sapply(pn,function(x) str_extract(x,"[[:digit:]]+"))
   wo=which(kegg_organism[,2]==input$Organism)
   organism=kegg_organism[wo,1] # Abbreviation of organism chosen
   code_path=paste0(organism,code_path)
  
   all_masses=round(as.numeric(raw[,2]),digits=4)
   row_list_non_peptide=setdiff(1:nrow(raw),row_list_peptide)# row in additional_data that not correspond to peptides
   masses_non_peptide=all_masses[row_list_non_peptide]
   masses_peptide=all_masses[row_list_peptide]
 
   list(code_path=code_path,pathway_names=pathways,mnp=masses_non_peptide,mp=masses_peptide)
 })
 

output$Kegg_pathway <- renderUI({
  selectInput('pathway', 'Choose your pathway:',load_KEGG()$pathway_names)
  #  code_path=paste0("sce",code_path)
})


output$image<-renderUI({
  
  i=which(load_KEGG()$pathway_names==input$pathway)
  code_path=load_KEGG()$code_path
  url=NULL
  
  ll=try(keggGet(code_path[i]),silent=T)
    if (class(ll)!="try-error"){
      cpds=try(names(ll[[1]]$COMPOUND),silent=T)
      if ((class(cpds)!="NULL") && (class(cpds)!="try-error")) {
        cpds_mass_list=sapply(cpds,function(x) as.numeric(keggGet(x)[[1]]$EXACT_MASS))
        cpds_mass_list=round(unlist(cpds_mass_list),4) # Round to 4 values after decimal
        mc=match(cpds_mass_list,load_KEGG()$mnp)
        mp=match(cpds_mass_list,load_KEGG()$mp)
        nodes1=names(cpds_mass_list)[which(!is.na(mp))]  # Peptides
        nodes2=names(cpds_mass_list)[which(!is.na(mc))] # No peptides but matched
        col1=rep('red',length(nodes1))
        col2=rep('yellow',length(nodes2))
        url= color.pathway.by.objects(code_path[i],c(nodes1,nodes2), c(col1,col2), c(col1,col2))
      }
    }
  tags$img(src=url)
  })
  
 
})
