enrich <-
function(id,nome_ris,percorso,file_set,nome_GO) {

require("csbl.go")
set.prob.table(type="enrichment",filename=file_set)
cat("set.prob.table DONE\n It has to be done in a new R session\n")
#leggo GO e tolgo id che non hanno GO

a=readLines(nome_GO)
a_sp=strsplit(a," ")
ff=file(paste("/tmp/GO_pul_",nome_ris,sep=""),"w")

for(i in 1:length(a)) {
  wh=which(id==a_sp[[i]][1])
  if(length(wh)!=0) {
    if(length(a_sp[[i]])!=1)
      cat(a[i],"\n",sep="",file=ff)
  }
}
close(ff)


ent=entities.from.text(paste("/tmp/GO_pul_",nome_ris,sep=""))
freq=frequencies(ent,p.values=T,fdr=T,normalize.one = TRUE)

if(dim(freq$BP)[1]!=0) {
  freq$BP[,"goid"]=ids.to.descs(freq$BP[,"goid"])
  write.table(freq$BP,quote=F,file=paste(percorso,nome_ris,"_freq_GO_BP",sep=""),sep="\t",row.names=TRUE)
  sel1=which(freq$BP[,"p.value"]<0.05)
  sel2=which(freq$BP[,"q.value"]<0.05)
  sel=intersect(sel1,sel2)
  tab=freq$BP[sel,c("proportion","priori")]
  go_tot=length(grep("GO:0008150",a))
  for(i in 1:dim(tab)[1]) {
    tab[i,2]=length(grep(rownames(tab)[i],a))/go_tot
  }
  #tolgo go con ref=0 perche' li mette ma non sono presenti in go slim
  rem=which(tab[,2]==0)
  #tolgo go che sono maggiori in ref rispetto a query
  rem2=which((tab[,1]-tab[,2])<0)
  if(length(rem2)!=0)
    rem=c(rem,rem2)
  if(length(rem)!=0)
    tab=tab[-rem,]
  colnames(tab)=c("query","ref")
  sel=rownames(tab)
  tab_fin=round(tab,3)
  desc=ids.to.descs(rownames(tab_fin))
  id_sel=readLines(paste("/tmp/GO_pul_",nome_ris,sep=""))
  ids=NULL
  for(i in 1:length(sel)) {
    ids=c(ids,paste(gsub(" GO:[G,O,0-9,:, ]+","",id_sel[grep(sel[i],id_sel)]),collapse=" "))
  }
  p_value=signif(freq$BP[rownames(tab_fin),"p.value"],3)
  q_value=signif(freq$BP[rownames(tab_fin),"q.value"],3)
  tab_fin=cbind(tab_fin,desc,p_value,q_value,ids)
  write.table(tab_fin,quote=F,file=paste(percorso,nome_ris,"_sel_GO_BP",sep=""),sep="\t",row.names=TRUE)
  tab=as.matrix(t(tab))
  barplot(tab,beside=TRUE,col=c("red","lightblue"),las=2, cex.names=0.8)
  savePlot(filename=paste(percorso,nome_ris,"_bartab_BP.jpg",sep=""),type="jpeg")
  graphics.off()
}

if(dim(freq$CC)[1]!=0) {
  freq$CC[,"goid"]=ids.to.descs(freq$CC[,"goid"])
  write.table(freq$CC,quote=F,file=paste(percorso,nome_ris,"_freq_GO_CC",sep=""),sep="\t",row.names=FALSE)
  sel1=which(freq$CC[,"p.value"]<0.05)
  sel2=which(freq$CC[,"q.value"]<0.05)
  sel=intersect(sel1,sel2)
  tab=freq$CC[sel,c("proportion","priori")]
  go_tot=length(grep("GO:0005575",a))
  for(i in 1:dim(tab)[1]) {
    tab[i,2]=length(grep(rownames(tab)[i],a))/go_tot
  }
  #tolgo go con ref=0 perche' li mette ma non sono presenti in go slim
  rem=which(tab[,2]==0)
  #tolgo go che sono maggiori in ref rispetto a query
  rem2=which((tab[,1]-tab[,2])<0)
  if(length(rem2)!=0)
    rem=c(rem,rem2)
  if(length(rem)!=0)
    tab=tab[-rem,]
  colnames(tab)=c("query","ref")
  sel=rownames(tab)
  tab_fin=round(tab,3)
  desc=ids.to.descs(rownames(tab_fin))
  id_sel=readLines(paste("/tmp/GO_pul_",nome_ris,sep=""))
  ids=NULL
  for(i in 1:length(sel)) {
    ids=c(ids,paste(gsub(" GO:[G,O,0-9,:, ]+","",id_sel[grep(sel[i],id_sel)]),collapse=" "))
  }
  p_value=signif(freq$CC[rownames(tab_fin),"p.value"],3)
  q_value=signif(freq$CC[rownames(tab_fin),"q.value"],3)
  tab_fin=cbind(tab_fin,desc,p_value,q_value,ids)
  write.table(tab_fin,quote=F,file=paste(percorso,nome_ris,"_sel_GO_CC",sep=""),sep="\t",row.names=TRUE)
  tab=as.matrix(t(tab))
  barplot(tab,beside=TRUE,col=c("red","lightblue"),las=2, cex.names=0.8)
  savePlot(filename=paste(percorso,nome_ris,"_bartab_CC.jpg",sep=""),type="jpeg")
  graphics.off()
}

if(dim(freq$MF)[1]!=0) {
  freq$MF[,"goid"]=ids.to.descs(freq$MF[,"goid"])
  write.table(freq$MF,quote=F,file=paste(percorso,nome_ris,"_freq_GO_MF",sep=""),sep="\t",row.names=FALSE)
  sel1=which(freq$MF[,"p.value"]<0.05)
  sel2=which(freq$MF[,"q.value"]<0.05)
  sel=intersect(sel1,sel2)
  tab=freq$MF[sel,c("proportion","priori")]
  go_tot=length(grep("GO:0003674",a))
  for(i in 1:dim(tab)[1]) {
    tab[i,2]=length(grep(rownames(tab)[i],a))/go_tot
  }
  #tolgo go con ref=0 perche' li mette ma non sono presenti in go slim
  rem=which(tab[,2]==0)
  #tolgo go che sono maggiori in ref rispetto a query
  rem2=which((tab[,1]-tab[,2])<0)
  if(length(rem2)!=0)
    rem=c(rem,rem2)
  if(length(rem)!=0)
    tab=tab[-rem,]
  colnames(tab)=c("query","ref")
  sel=rownames(tab)
  tab_fin=round(tab,3)
  desc=ids.to.descs(rownames(tab_fin))
  id_sel=readLines(paste("/tmp/GO_pul_",nome_ris,sep=""))
  ids=NULL
  for(i in 1:length(sel)) {
    ids=c(ids,paste(gsub(" GO:[G,O,0-9,:, ]+","",id_sel[grep(sel[i],id_sel)]),collapse=" "))
  }
  p_value=signif(freq$MF[rownames(tab_fin),"p.value"],3)
  q_value=signif(freq$MF[rownames(tab_fin),"q.value"],3)
  tab_fin=cbind(tab_fin,desc,p_value,q_value,ids)
  write.table(tab_fin,quote=F,file=paste(percorso,nome_ris,"_sel_GO_MF",sep=""),sep="\t",row.names=TRUE)
  tab=as.matrix(t(tab))
  barplot(tab,beside=TRUE,col=c("red","lightblue"),las=2, cex.names=0.8)
  savePlot(filename=paste(percorso,nome_ris,"_bartab_MF.jpg",sep=""),type="jpeg")
  graphics.off()
}

}
