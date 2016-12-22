options(warn=-1)


##########################
#### Auxillary functions:
###########################

cut_vector<-function(ll){
  # ll is a numeric vector
  output=list()
  if (length(ll)>1){
  k=1
  output[[k]]=ll[1]
  for (i in 2:length(ll)){
    if (ll[i]-ll[i-1]==1){
      output[[k]]=c(output[[k]],ll[i])}
    else {k=k+1
    output[[k]]=ll[i]}}}
  else{output[[1]]=ll}
  return(output)
}
  
composition_formula<-function(formula,elements){
  nb_elements=rep(0,length(elements))
  splitted=strsplit(formula,'')[[1]]
  
  is_letter=which(!grepl('\\d',splitted))
  is_letter=cut_vector(is_letter)
  
  is_number=which(grepl('\\d',splitted))
  is_number=cut_vector(is_number)
  
  for (i in 1:length(is_letter)){
    w=which(elements==paste0(splitted[is_letter[[i]]],collapse=""))
    nb_elements[w]=as.numeric(paste0(splitted[is_number[[i]]],collapse=""))}
  return(nb_elements)}

filter_graph<-function(g,from_list,to_list){ # If PLA->PL, PLA->L and PL->L coexists, filter PLA->L)
  valid=c()
  for (i in 1:length(to_list)){
    l=shortest_paths(g,from_list[i],to_list[i])
    g1=delete_edges(g,paste0(from_list[i],"|",to_list[i])) # Remove edges
    l2=shortest_paths(g1,from_list[i],to_list[i])
    if (length(l2$vpath[[1]])==0){valid=c(valid,i)} # No more edges beween
  }
  return(valid)}

# Check if the annotation is an amino acid
check_free<-function(annotation,elements){
  tmp=composition_formula(annotation,elements)
  return(sum(tmp)!=1)
} # True only when its a peptide

# Check if both elements of l are in the node list
check_in<-function(l,nodes){
  return((l[1]%in%nodes)&(l[2]%in%nodes))}

# The following function request the job from decomp server

send_curl<-function(mass_list,monomers,tol,DecompID,dplace){

# Prepare Json style input & parameters:  

if (DecompID=="Job ID"){
  
  header="{ \"decomp_input_real_alphabet\":\"# monomer masses, monoisotopic distribution"
  
  units=paste0("\\r\\n",monomers[,1]," ",round(monomers[,2],dplace),collapse="")
  
  middle="\", \"decomp_input_real_masses\":\""

  masses=paste0(mass_list[1:(length(mass_list)-1)],"\\r\\n",collapse="")
  
  masses=paste0(masses,mass_list[length(mass_list)],collapse="")
  
  ender="\", \"paramset\":{ \"decomp_masses_masserrorunit\":\"Da\", \"decomp_masses_massdistribution\":\"mono\", \"decomp_filtering_deviation\":\"true\", \"decomp_filtering_best\":\"20\", \"decomp_filtering_chemicallyplausible\":\"false\", \"decomp_masses_allowedmasserror\":\""
  
  ender2=paste0(tol,"\"} }")
  
  input_str=paste0(header,units,middle,masses,ender,ender2,collapse="")

# 1: Post data

#curl -X POST -d @[[INPUT]] http://bibiserv2.cebitec.uni-bielefeld.de:80/rest/decomp/decomp_decompose_reals/request -H "Content-Type: application/json"

  b1="http://bibiserv2.cebitec.uni-bielefeld.de:80/rest/decomp/decomp_decompose_reals/request"
  p1=postForm(b1,.opts=list(postfields=input_str, httpheader="Content-Type: application/json",useragent = "RCurl"))

# 2: Check status
  b2="http://bibiserv2.cebitec.uni-bielefeld.de:80/rest/decomp/decomp_decompose_reals/statuscode"
  id=p1[1]
  ptm0 <- proc.time()
  p2="0"
  time_passed=0
  while (p2!="600" & time_passed<500){
    p2=postForm(b2,.opts=list(postfields=id, httpheader="Content-Type: text/plain",useragent = "RCurl"))
    ptm1 <- proc.time()-ptm0  
    time_passed=ptm1[3]
  }
  
# 3: Output

#curl -X POST -d [[ID]] http://bibiserv2.cebitec.uni-bielefeld.de:80/rest/decomp/decomp_decompose_reals/response -H "Content-Type: text/plain";echo
  if (p2=="600"){ # If the job is finished with success
    b3="http://bibiserv2.cebitec.uni-bielefeld.de:80/rest/decomp/decomp_decompose_reals/response"
    p3=postForm(b3,.opts=list(postfields=id, httpheader="Content-Type: text/plain",useragent = "RCurl"))
  }
  if (p2!="600") {p3="No response from the server"}} # Send error message if Decomp faied to find the results
 
else{
  b3="http://bibiserv2.cebitec.uni-bielefeld.de:80/rest/decomp/decomp_decompose_reals/response"
  id=DecompID
  p3=postForm(b3,.opts=list(postfields=DecompID, httpheader="Content-Type: text/plain",useragent = "RCurl"))
  }
  return(list(p3=p3,id=id))
}

### Create KEGG links

createLink <- function(val) {
  link_list=c()
  for (i in 1:length(val)){
    if (val[i]=="0"){link_list=c(link_list,	'NA')}
    else {
      cpds=strsplit(val[i], " , ")[[1]] # All cpds
      cpds=paste0(cpds,collapse="+")
      site=paste0('<a href="http://www.genome.jp/dbget-bin/www_bget?',cpds,'">',cpds,'</a>')
    link_list=c(link_list,site)}}
   return(link_list) 
  }
 #'<a href="http://rstudio.com">RStudio</a>'

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

decimalnumcount<-function(x){stopifnot(class(x)=="character")
  x<-gsub("(.*)(\\.)|([0]*$)","",x)
  nchar(x)
}

