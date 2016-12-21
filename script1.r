
options(warn=-1)

####################
#### Main functions:
###################

# Script for data annotation:

peptide_annotation<-function(raw_data,additional_data,results,tol){
  
  peptide_annotated=c() # List of all peptide annotations (including no-annotation & duplex annotations)
  nb_peptide_annotated=c() # List of number of peptide annotations for each mass
  
  unique_peptide_annotated=c() # List of valid peptides with unique annotation
  unique_row=c() # Selected row that represent an unique annotation of peptides
  
  # Split results to different masses
  results0= results
  results=strsplit(results0,split="\n# done\n")
  results=results[[1]][1]
  results=strsplit(results,split="\n# mass")
  results=results[[1]][3:length(results[[1]])]
  
  for (i in 1:length(results)){
    
    annotations=strsplit(results[i],"\n")
    annotations=annotations[[1]]
    if (length(annotations)>1){
      
      possible_list=c() # All possible annotations for the same mass
      for (k in 2:length(annotations)){
        error=as.numeric(gsub(".*\\((.*)\\).*", "\\1",  annotations[k]))
        
        if (abs(error)< tol){ # Computational error filtering
          letters=strsplit(annotations[k],split="\\(")[[1]][1]
          possible_list=c(possible_list,gsub(" ","",letters))}}
      
      if (length(possible_list)==0){ # if no peptide annoation possible
        nb_peptide_annotated=c(nb_peptide_annotated,0)
        peptide_annotated=c(peptide_annotated,"No Peptide Found")
      }
      
      else if (length(possible_list)==1){ # If only one peptide annotation possible
        
        peptide_annotated=c(peptide_annotated,possible_list)
        nb_peptide_annotated=c(nb_peptide_annotated,1)
        unique_peptide_annotated=c(unique_peptide_annotated,possible_list) # List of valid peptides with unique annotation
        unique_row=c(unique_row,i)}
      
      else {peptide_annotated=c(peptide_annotated,paste(possible_list,collapse=";"))
      nb_peptide_annotated=c(nb_peptide_annotated,length(possible_list))}
    }
    else{
      peptide_annotated=c(peptide_annotated,"No Peptide Found")
      nb_peptide_annotated=c(nb_peptide_annotated,0)
    }
  }
      # If multiple possible annotations

 raw_data_additional=cbind(additional_data,NBP=nb_peptide_annotated,Peptide=peptide_annotated)
 unique_data_annotated=cbind(raw_data[unique_row,],Peptide=unique_peptide_annotated)
 unique_data_additional=cbind(additional_data[unique_row,],Peptide=unique_peptide_annotated)

 return(list(all=raw_data_additional,unique=unique_data_annotated,unique_add=unique_data_additional))}

# Script for edge & node calculation

network_generator<-function(unique,unique_add,cor_min,index,elements){
  
  I_matrix=unique[,3:(ncol(unique)-1)]
  class(I_matrix)='numeric'
  
  elements=sort(elements)
  ID_list=as.numeric(unique[,1])
  from_list=c()
  to_list=c()
  diff_list=c()
  cor_list=c()
  
  peptides=unique[,ncol(unique)]
  peptide_list=lapply(peptides,composition_formula,elements)

  for (p in 1:(length(peptide_list)-1)){
    for (q in (p+1):length(peptide_list)){
      dis=peptide_list[[q]]-peptide_list[[p]]
      coef=cor(I_matrix[p,],I_matrix[q,],method="spearman")
     if (is.na(coef)){coef=-1}
      if (all(dis>=0) && coef>=cor_min){
        from_list=c(from_list,ID_list[q])
        to_list=c(to_list,ID_list[p])
        diff_names=elements[which(dis>0)]
        diff_list=c(diff_list,paste(paste0(diff_names,dis[which(dis>0)]),collapse=''))
        cor_list=c(cor_list,coef)
     }
   }
  }
  
  vectorized=as.vector(rbind(from_list,to_list))
  g=graph(vectorized)
  cor_list=as.numeric(cor_list)
  
  valid=filter_graph(g,from_list,to_list) # Filter short chains 
  from_list=from_list[valid]
  to_list=to_list[valid]
  diff_list=diff_list[valid]
  cor_list=cor_list[valid]
  
  vectorized=as.vector(rbind(from_list,to_list))
  g=graph(vectorized) # igraph object
  
  links=data.frame(from = from_list, to = to_list)
  links$arrows <- "to"
  links$smooth <- FALSE    # should the edges be curved?
  links$shadow <- FALSE    # edge shadow
  links$title=diff_list
  
  unique_ID=sort(unique(c(from_list,to_list)))
  ID_list=as.numeric(unique_add[,1])
  matched_row=match(unique_ID,ID_list)
  unique_add2=unique_add[matched_row,]
  
  nodes <- data.frame(id = as.numeric(unique_add2[,1]))
  nodes$shape='dot'
  nodes$shadow <- TRUE 
  nodes$title <- unique_add2[,index]
  nodes$label <- unique_add2[,ncol(unique_add2)]
  nodes$color.background='blue'
  
  edges=cbind(from_list,to_list,diff_list)
  colnames(edges)=c("Source","Target","Loss")
  
return(list(links=links,nodes=nodes,from=from_list,to=to_list,diff_list=diff_list,edges=edges,cor=cor_list,g=g))}

# Find common pattern modules

common_pattern_single<-function(network,node_id,mode,orders){
  
  g=network$g
  
#  L0=0
#  L1=1
#  d=1
#  nodes_list=node_id
  
#  while (L1>L0){
#    L0=length(nodes_list)
#    connected_nodes=ego(g,d,mode=mode,nodes=node_id,mindist=0)
#    nodes_list=connected_nodes[[1]]
#    L1=length(nodes_list)
#    d=d+1}

  connected_nodes=ego(g,orders,mode=mode,nodes=node_id,mindist=0)
  nodes_list=connected_nodes[[1]]
  valid_links=which(apply(rbind(network$from,network$to),2,check_in,nodes_list))
  
  new_links=network$links[valid_links,]
  from_list=network$from[valid_links]
  to_list=network$to[valid_links]
  diff_list=network$diff_list[valid_links]
  cor_list=network$cor[valid_links]
  
  vectorized=as.vector(rbind(from_list,to_list))
  g=graph(vectorized) # igraph object
  
  unique_ID=sort(unique(c(from_list,to_list)))
  ID_list=as.numeric(network$nodes[,1])
  matched_row=match(unique_ID,ID_list)
  new_nodes <- network$nodes[matched_row,]

  edges=cbind(from_list,to_list,diff_list)
  colnames(edges)=c("Source","Target","Loss")
  
  # edge shadow
  
  return(list(links=new_links,nodes=new_nodes,from=from_list,to=to_list,diff_list=diff_list,edges=edges,cor=cor_list,g=g))
}

# Find all degradation chains

degrade_chain_all<-function(network,degree_size){
  
  g=network$g
  nodes=network$nodes
  start_labels=c()
  paths=list()
  d=1
  for (n in 1:(nrow(nodes)-1)){
    for (m in ((n+1):nrow(nodes))){
      checked_nodes=all_simple_paths(g,from=nodes[m,1],to=nodes[n,1])
      if (length(checked_nodes)>=1){
        chain_lengths=lapply(checked_nodes,length) # Length of all found connections between nodes
        wm=which(chain_lengths>degree_size)
        if (length(wm)>0){
          for (i in 1:length(wm)){ 
            chain_nodes=checked_nodes[[wm[i]]] # One of the long chains
            chain_nodes=sort(chain_nodes,decreasing=T)
            paths[[d]]=chain_nodes
            valid=match(chain_nodes[1],nodes$id)
            start_labels=c(start_labels,paste0(nodes$label[valid],"-Chain-",d,"-Length-",length(chain_nodes)-1))
            d=d+1}
      }}
    }
  }
 valid=which(!duplicated(paths)) # no duplicated paths
 start_labels=start_labels[valid] 
 paths_new=list()
 for (i in 1:length(valid)){paths_new[[i]]=paths[[valid[i]]]}
 return(list(paths=paths_new,start_labels=start_labels)) 
}

generate_degrade_chain<-function(network,path){
  
  valid_links=which(apply(rbind(network$from,network$to),2,check_in,path))
  new_links=network$links[valid_links,]
  from_list=network$from[valid_links]
  to_list=network$to[valid_links]
  cor_list=network$cor[valid_links]
  diff_list=network$diff_list[valid_links]
  
  vectorized=as.vector(rbind(from_list,to_list))
  g=graph(vectorized) # igraph object
  
  unique_ID=sort(unique(c(from_list,to_list)))
  ID_list=as.numeric(network$nodes[,1])
  matched_row=match(unique_ID,ID_list)
  new_nodes <- network$nodes[matched_row,]
  
  edges=cbind(from_list,to_list,diff_list)
  colnames(edges)=c("Source","Target","Loss")
  
  # edge shadow
  
  return(list(links=new_links,nodes=new_nodes,from=from_list,to=to_list,diff_list=diff_list,edges=edges,cor=cor_list,g=g))
}

# Remove amino acids from network
filter_amino_acid<-function(network,elements){ 
  annotation_list=network$nodes$label
  wp=which(sapply(annotation_list,check_free,elements))
  new_nodes <- network$nodes[wp,]
  valid_links=which(apply(rbind(network$from,network$to),2,check_in,new_nodes$id))
  
  new_links=network$links[valid_links,]
  from_list=network$from[valid_links]
  to_list=network$to[valid_links]
  cor_list=network$cor[valid_links]
  diff_list=network$diff_list[valid_links]
  
  vectorized=as.vector(rbind(from_list,to_list))
  g=graph(vectorized) # igraph object
  
  edges=cbind(from_list,to_list,diff_list)
  colnames(edges)=c("Source","Target","Loss")
  
return(list(links=new_links,nodes=new_nodes,from=from_list,to=to_list,diff_list=diff_list,edges=edges,cor=cor_list,g=g))}

