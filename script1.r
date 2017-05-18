
options(warn=-1)

##############################
#Scripts for mass decomposition
##############################

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


# Script for data annotation from DECOMP results:

peptide_annotation<-function(raw_data,additional_data,results,tol,monomers,condensation){
  
  peptide_annotated=c() # List of all peptide annotations (including no-annotation & duplex annotations)
  nb_peptide_annotated=c() # List of number of peptide annotations for each mass
  error_annotated=c() # List of ppm errors for each mass
  
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
          possible_list=c(possible_list,gsub(" ","",letters))}
      }
      
      if (length(possible_list)==0){ # if no peptide annoation possible
        nb_peptide_annotated=c(nb_peptide_annotated,0)
        peptide_annotated=c(peptide_annotated,"No Peptide Found")
        error_annotated=c(error_annotated,"NA")
      }
      
      else if (length(possible_list)==1){ # If only one peptide annotation possible
        
        peptide_annotated=c(peptide_annotated,possible_list)
        nb_peptide_annotated=c(nb_peptide_annotated,1)
        unique_peptide_annotated=c(unique_peptide_annotated,possible_list) # List of valid peptides with unique annotation
        unique_row=c(unique_row,i)
        ppm_error=ppm_calc(possible_list,monomers,condensation,raw_data[i,2])
        error_annotated=c(error_annotated,ppm_error)}
      
      else { # If multiple annotation possible
        peptide_annotated=c(peptide_annotated,paste(possible_list,collapse=";"))
        nb_peptide_annotated=c(nb_peptide_annotated,length(possible_list))
        error_list=c()
        for (t in 1:length(possible_list)){
          ppm_error=ppm_calc(possible_list[t],monomers,condensation,raw_data[i,2])
          error_list=c(error_list,ppm_error)}
        ppm_error=paste0(error_list,collapse=";")
        error_annotated=c(error_annotated,ppm_error)}
    }
    else{
      peptide_annotated=c(peptide_annotated,"No Peptide Found")
      nb_peptide_annotated=c(nb_peptide_annotated,0)
      error_annotated=c(error_annotated,"NA")
    }
  }
      # If multiple possible annotations
 raw_data=cbind(raw_data,NBP=nb_peptide_annotated,Peptide=peptide_annotated,PPM=error_annotated)
 raw_data_additional=cbind(additional_data,NBP=nb_peptide_annotated,Peptide=peptide_annotated,PPM=error_annotated)
 unique_data_annotated=raw_data[unique_row,]
 unique_data_additional=cbind(additional_data[unique_row,],Peptide=unique_peptide_annotated,PPM=error_annotated[unique_row])

 return(list(all=raw_data,all_add=raw_data_additional,unique=unique_data_annotated,unique_add=unique_data_additional))}

# Script for data annotation from recursive functinon:


peptide_annotation_slow<-function(masslist,dplace,monomers,raw_data,additional_data,tol,condensation){
  
  peptide_annotated=c() # List of all peptide annotations (including no-annotation & duplex annotations)
  nb_peptide_annotated=c() # List of number of peptide annotations for each mass
  unique_row=c() # Selected row that represent an unique annotation of peptides
  error_annotated=c()
  
  monomers[,2]=round(monomers[,2],dplace)
  
  if (tol<1e-4) {appro_dplace=3}
  if ((tol>=1e-4) & (tol<=1e-1)) {appro_dplace=-round(log10(tol))}
  if (tol>1e-1) {appro_dplace=0}

  masslist_arrondi=round(masslist,appro_dplace)*(10^appro_dplace)
  monomers_arrondi=round(monomers[,2],appro_dplace)*(10^appro_dplace) # Make every value integer
  
  for (i in 1:length(masslist_arrondi)){
  #  print(masslist_arrondi[i])
    mtry=c(masslist_arrondi[i]-1,masslist_arrondi[i],masslist_arrondi[i]+1) # Test different for rounding effect
    possible_combinations=c()
    for (mass in mtry){
      tmp=capture.output(pay(mass,monomers_arrondi,c()))
      if (length(tmp)>0){
        for (k in 1:length(tmp)){
          kt=strsplit(tmp[k]," ")[[1]]
          kt=kt[3:length(kt)]
          possible_combinations=rbind(possible_combinations,as.numeric(kt))}}
      if (!is.null(possible_combinations)){ # Decomposition succeeded
       possible_combinations=unique(possible_combinations)
       calculated_masslist=possible_combinations%*%monomers[,2]
       valid=which(abs(t(calculated_masslist)-masslist[i])<=tol)# Confirm correct decompositions
       possible_combinations=matrix(possible_combinations[valid,],ncol=nrow(monomers))
    if (nrow(possible_combinations)>10) break}}   # Until 10 peptides only
    
    if (!is.null(possible_combinations)){
      if (nrow(possible_combinations)>0){
      possible_peptides=c()
      error_list=c()
      for (r in 1:nrow(possible_combinations)){
        position=which(possible_combinations[r,]>0)
        peptide=paste(paste0(monomers[position,1],possible_combinations[r,position]),collapse="")
        ppm_error=ppm_calc(peptide,monomers,condensation,raw_data[i,2])
        possible_peptides=c(possible_peptides,peptide)
        error_list=c(error_list,ppm_error)}
      
      possible_peptides=paste(possible_peptides,collapse=" ; ")
      error_list=paste(error_list,collapse=" ; ")
      peptide_annotated=c(peptide_annotated,possible_peptides)
      error_annotated=c(error_annotated,error_list)
      nb_peptide_annotated=c(nb_peptide_annotated,length(valid))}
    
    if (nrow(possible_combinations)==1) {unique_row=c(unique_row,i)}
     
    if (nrow(possible_combinations)==0){   
        peptide_annotated=c(peptide_annotated,"No Peptide Found")
        nb_peptide_annotated=c(nb_peptide_annotated,0)
        error_annotated=c(error_annotated,"NA")}}
    
    else{peptide_annotated=c(peptide_annotated,"No Peptide Found")
    nb_peptide_annotated=c(nb_peptide_annotated,0)
    error_annotated=c(error_annotated,"NA")}}

   raw_data=cbind(raw_data,NBP=nb_peptide_annotated,Peptide=peptide_annotated,PPM=error_annotated)
   raw_data_additional=cbind(additional_data,NBP=nb_peptide_annotated,Peptide=peptide_annotated,PPM=error_annotated)
   unique_data_annotated=raw_data[unique_row,]
   unique_data_additional=cbind(additional_data[unique_row,],Peptide=peptide_annotated[unique_row],PPM=error_annotated[unique_row])
  
   return(list(all=raw_data,all_add=raw_data_additional,unique=unique_data_annotated,unique_add=unique_data_additional))}

# Recursive function:
pay<-function(Amount,billList,billCount){
  billList=sort(billList,decreasing=T)       
  maxN=floor(Amount/billList[1])
  if (length(billList)==1){
    reminder=Amount%%billList
    if (reminder==0){
      count=c(billCount,maxN)
      print(count)
    }
  }
  if (length(billList)>1){
    for (i in seq(maxN,0,by=-1)){
      reminder=Amount-i*billList[1]
      tmp=pay(reminder,billList[2:length(billList)],c(billCount,i))
      }
  }
}

annotate_example<-function(url1,url2,monomers,tol2){
  
  raw_data=data.matrix(read.table(url1,sep='\t',dec=',',header=T,stringsAsFactors = F))
  
  additional_data=data.matrix(read.table(url2,sep='\t',dec=',',header=T,stringsAsFactors = F))
  
  dplace=decimalnumcount(as.character(raw_data[1,2])) # decimal place
  
  mass_list=raw_data[,2]-round(monomers$condensation,dplace)
  
  tol1=0.01 # Default tolerance for Decomp (Theoritical)
  
  results=send_curl(mass_list,monomers$tab,tol1,"Job ID",dplace) # bielefeld results
  
  output_message="Decomposition succeeded! please download decomposition results!"
  output_massage=c(output_message,paste0("You Job id on DECOMP server is: ",results$id))
  annotated=peptide_annotation(raw_data,additional_data,results$p3,tol2)
  
  return(list(annotated=annotated,output_message=output_message,raw_data=raw_data,dplace=dplace,tol2=tol2))
}
















