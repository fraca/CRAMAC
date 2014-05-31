redGO <-
function(file_in,nome,file_obo,onto="") {

a=readLines(file_in)

a=strsplit(a," ")

b=readLines(file_obo)

names(b)=substr(b,1,10)

b=gsub("^[GO,:,0-9,=,>, ]+ // ","",b)
#b=strsplit(b," // ")


ris=NULL
cat("It take some time (",length(a),") wait...\n")

for(i in 1:length(a)) {
  gon=unlist(strsplit(b[a[[i]][2]]," "))
  if(length(gon)==0)
    stop("AHAHAHA problema")
  for(j in 1:length(gon)) {
    if(onto!="") {
      if(ids.to.ontologies(gon[j])==onto)
	ris=c(ris,paste(a[[i]][1],gon[j]))
    } else {
      ris=c(ris,paste(a[[i]][1],gon[j]))
    }
  }
  #cat(i,length(a),"\n")
}
ris=unique(ris)
cat(ris,sep="\n",file=nome)


}
