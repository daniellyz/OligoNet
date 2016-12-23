source("script1.r")
source("script2.r")

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
