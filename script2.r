options(warn=-1)

###############################
# Script for network construction
###############################

# Script for edge & node calculation

network_generator<-function(unique,unique_add,cor_min,cor_max,elements,delete_triangle,isolated,filter_aa){
  
  ranked=order(unique[,2])
  unique=unique[ranked,]
  unique_add=unique_add[ranked,]
  
  I_matrix=unique[,3:(ncol(unique)-3)]
 # class(I_matrix)='numeric'
  I_matrix=data.matrix(I_matrix)  
  
  elements=sort(elements)
  ID_list=as.character(unique[,1])
  ID_list0=ID_list
  
  from_list=c()
  to_list=c()
  diff_list=c()
  cor_list=c()
  cor_invalid_list=c() # Correlation calculation impossible

  peptides=unique[,"Peptide"]

  if (filter_aa==1){
    is_peptide_list=which(sapply(peptides,check_free,elements))
    peptides=peptides[is_peptide_list]
    ID_list=ID_list[is_peptide_list]
    I_matrix=I_matrix[is_peptide_list,]
    unique=unique[is_peptide_list,]
    unique_add=unique_add[is_peptide_list,]}
  
  peptide_list=lapply(peptides,composition_formula,elements)

  for (p in 1:(length(peptide_list)-1)){
    for (q in (p+1):length(peptide_list)){
      dis=peptide_list[[q]]-peptide_list[[p]]
      coef=cor(I_matrix[p,],I_matrix[q,],method="spearman",use="pairwise.complete.obs")
      if (is.na(coef)){coef=-100} # no valid correlations
      if (all(dis>=0) && coef>=cor_min && coef<=cor_max){
        from_list=c(from_list,ID_list[q])
        to_list=c(to_list,ID_list[p])
        diff_names=elements[which(dis>0)]
        if (all(dis==0)){diff_list=c(diff_list,"Isomers")}
        else{diff_list=c(diff_list,paste(paste0(diff_names,dis[which(dis>0)]),collapse=''))}
        cor_list=c(cor_list,coef)
        if (is.na(coef)){cor_invalid_list=c(cor_invalid_list,1)} # no valid correlations
        else{cor_invalid_list=c(cor_invalid_list,0)}
      }
    }
  }
  
  vectorized=as.vector(rbind(from_list,to_list))
  g=graph(vectorized) # igraph object
  cor_list=as.numeric(cor_list)
  
  if (delete_triangle==1){
    valid=filter_graph(g,from_list,to_list) # Filter short chains 
    from_list=from_list[valid]
    to_list=to_list[valid]
    diff_list=diff_list[valid]
    cor_list=cor_list[valid]}
  
  vectorized=as.vector(rbind(from_list,to_list))
  g=graph(vectorized) # igraph object
  
  links=data.frame(from = from_list, to = to_list,stringsAsFactors = F)
  links$arrows <- "to"
  links$smooth <- FALSE    # should the edges be curved?
  links$shadow <- FALSE    # edge shadow
  links$title=diff_list
  
  network_ID=unique(c(from_list,to_list)) # Only ID connected
  matched_row=match(network_ID,ID_list)
  
  if (isolated==0){unique2=unique[matched_row,]} # only the nodes in the network
  else {unique2=unique} # All including isolated nodes
  
  nodes <- data.frame(id = as.character(unique2[,1]),stringsAsFactors = F)
  nodes$shape='dot'
  nodes$shadow<-TRUE 
  nodes$label <- unique2[,"Peptide"]
  
  from_valid=match(from_list,ID_list) # Which row in unique
  to_valid=match(to_list,ID_list)
  id_edges=paste0("E",1:length(from_list))
  cor_list=round(cor_list,3)
  cor_list[which(cor_invalid_list==1)]="Not possible to calculate"
  edges=cbind(id_edges,from_list,unique[from_valid,"Peptide"],to_list,unique[to_valid,"Peptide"],paste0("- ",diff_list),cor_list)

  colnames(edges)=c("Edges ID","Source ID","Source Peptide","Target ID","Target Peptide","Amino acid(s) Loss","Statistical Correlations")
  
  matched=match(as.character(unique2[,1]),ID_list0) # matched used for labeling during network visualization

  return(list(links=links,nodes=nodes,from=from_list,to=to_list,diff_list=diff_list,edges=edges,cor=cor_list,g=g,matched=matched))}

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
  if (length(connected_nodes)==0){}
  else{
  nodes_list=connected_nodes[[1]]
  valid_links=which(apply(rbind(network$from,network$to),2,check_in,names(nodes_list)))
  
  new_links=network$links[valid_links,]
  from_list=network$from[valid_links]
  to_list=network$to[valid_links]
  diff_list=network$diff_list[valid_links]
  cor_list=network$cor[valid_links]
  new_edges=network$edges[valid_links,]
  
  vectorized=as.vector(rbind(from_list,to_list))
  g=graph(vectorized) # igraph object
  
  unique_ID=unique(c(from_list,to_list))
  ID_list=network$nodes[,1]
  matched_row=match(unique_ID,ID_list)
  new_nodes <- network$nodes[matched_row,]
  
  # edge shadow
  
  return(list(links=new_links,nodes=new_nodes,from=from_list,to=to_list,diff_list=diff_list,edges=new_edges,cor=cor_list,g=g))}
}

# Find all degradation chains

degrade_chain_all<-function(network,degree_size){
  
  g=network$g
  nodes=network$nodes
  start_labels=c()
  paths=list()
  paths_new=list()
  d=1
  
  for (n in 1:(nrow(nodes)-1)){
    for (m in ((n+1):nrow(nodes))){
      checked_nodes=try(all_simple_paths(g,from=nodes[m,1],to=nodes[n,1]),silent=T)
      if ((length(checked_nodes)>=1) & (class(checked_nodes)!="try-error")){
        chain_lengths=lapply(checked_nodes,length) # Length of all found connections between nodes
        wm=which(chain_lengths>degree_size)
        if (length(wm)>0){
          for (i in 1:length(wm)){ 
            chain_nodes=checked_nodes[[wm[i]]] # One of the long chains
            chain_nodes=sort(chain_nodes,decreasing=T)
            paths[[d]]=chain_nodes
            valid=match(names(chain_nodes[1]),nodes$id)
            start_labels=c(start_labels,paste0(nodes$label[valid],"-Chain-",d,"-Length-",length(chain_nodes)-1))
            d=d+1}
        }}
    }
  }
  
  if (length(paths)>0){
  valid=which(!duplicated(paths)) # no duplicated paths
  start_labels=start_labels[valid] 
  for (i in 1:length(valid)){paths_new[[i]]=paths[[valid[i]]]}}
  return(list(paths=paths_new,start_labels=start_labels)) 
}

generate_degrade_chain<-function(network,path){
  
  valid_links=which(apply(rbind(network$from,network$to),2,check_in,names(path)))
  new_links=network$links[valid_links,]
  from_list=network$from[valid_links]
  to_list=network$to[valid_links]
  cor_list=network$cor[valid_links]
  diff_list=network$diff_list[valid_links]
  new_edges=network$edges[valid_links,]
  
  vectorized=as.vector(rbind(from_list,to_list))
  g=graph(vectorized) # igraph object
  
  unique_ID=unique(c(from_list,to_list))
  ID_list=network$nodes[,1]
  matched_row=match(unique_ID,ID_list)
  new_nodes <- network$nodes[matched_row,]
  
  # edge shadow
  
  return(list(links=new_links,nodes=new_nodes,from=from_list,to=to_list,diff_list=diff_list,edges=new_edges,cor=cor_list,g=g))
}


