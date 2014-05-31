geni2coo <-
function(id,db,transcript=FALSE) {

coo_tot=paste(db,"_coo",sep="")

coo_tot=read.table(coo_tot,sep=";")
coo_tot[,1]=gsub(" seq=cdna","",coo_tot[,1])
coo_tot[,1]=gsub(">","",coo_tot[,1])
coo_tot[,2]=gsub(" coord=","",coo_tot[,2])
#coo_tot[,2]=gsub(":-?1$","",coo_tot[,2])
coo_tot[,3]=gsub(" parent_gene=","",coo_tot[,3])

if(transcript) {
 sel=NULL
 for(i in 1:length(id)) {
  se=which(coo_tot[,1]==id[i])
  if(length(se)!=0)
    sel=c(sel,se)
 }
} else {
 sel=NULL
 for(i in 1:length(id)) {
  se=which(coo_tot[,3]==id[i])
  if(length(se)!=0)
    sel=c(sel,se)
 }
}


tab_ris=coo_tot[sel,]
tab_ris=tab_ris[order(tab_ris[,2]),]

return(tab_ris)

}
