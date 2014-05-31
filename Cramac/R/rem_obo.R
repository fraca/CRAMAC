rem_obo <-
function(file_in,obs,file_out) {

agri=readLines(file_in)
obs=read.table(obs,sep="\t",stringsAsFactors=FALSE,header=TRUE)



#tolgo obosoleti o li cambio

rem=NULL
agg=NULL
for(i in 1:dim(obs)[1]) {
  aa=grep(obs[i,1],agri)
  if(length(aa)!=0) {
    rem=c(rem,aa)
    if(!is.na(obs[i,2])) {
      gon=unlist(strsplit(obs[i,2]," "))
      idn=gsub(" GO:[0-9]+","",agri[aa])
      cat(length(gon),length(idn),i,"\n")
      for(j in 1:length(gon)) {
	for(z in 1:length(idn)) {
	  agg=c(agg,paste(idn[z],gon[j],sep=" "))
	}
      }      
    }
  }
}
agri=agri[-rem]
agri=unique(c(agri,agg))

#a mano sistemo GO_term non presenti in cbsl.go sono tutti circa phosphatase activity
agri=gsub("GO:0038036","GO:0003376",agri) # da MF a BP
agri=gsub("GO:0052742","GO:0046854",agri) # da MF a BP
agri=gsub("GO:0052813","GO:0046854",agri) # da MF a BP
agri=gsub("GO:0052744","GO:0016791",agri)
agri=gsub("GO:0052745","GO:0016791",agri)
cat(agri,sep="\n",file=file_out)
}
