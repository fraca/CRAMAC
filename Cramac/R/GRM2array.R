GRM2array <-
function(id_maize,tab_conf,score_length=NULL) {

tab_c=read.table(tab_conf,header=TRUE,stringsAsFactors=FALSE)

if(!is.null(score_length))
  tab_c=tab_c[which(tab_c[,3]>score_length),]

tab_ris=tab_c[1,]
rimossi=NULL
for(i in 1:length(id_maize)) {
  nn=grep(id_maize[i],tab_c[,2])
  if(length(nn)!=0) {
    for(j in 1:length(nn)) {
      tab_ris=rbind(tab_ris,tab_c[nn[j],])
    }
  } else {
    rimossi=c(rimossi,id_maize[i])
  }
}

cat("Sono stati rimossi",length(rimossi),"su",length(id_maize),"id GRM che non sono presenti in tabella (/tmp/rimossi)\n")
cat(rimossi,sep="\n",file="/tmp/rimossi")
tab_ris=tab_ris[-1,]


return(tab_ris)

}
