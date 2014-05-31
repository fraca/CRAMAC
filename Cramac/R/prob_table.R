prob_table <-
function(nome,tab_in,dir) {

id=unique(tab_in[,1])
ff=file(paste(dir,nome,sep=""),"w")
for(i in 1:length(id)) {
   cat(id[i],tab_in[which(tab_in[,1]==id[i]),2],sep=" ",file=ff)
   cat("\n",file=ff)
 
 }
close(ff)

ent=entities.from.text(paste(dir,nome,sep=""))

protab=entities.to.prob.table(ent)
write.table(protab,file=paste(dir,"prob_tab_sim_",nome,sep=""),col.names=F,row.names=F)

protab_en=entities.to.enrichment.probs(ent)
write.table(protab_en,file=paste(dir,"prob_tab_enr_",nome,sep=""),col.names=F,row.names=F)

}
