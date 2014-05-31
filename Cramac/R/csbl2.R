csbl2 <-
function(nome_ris,percorso,n_lim=5,file_set,nome_GO) {


##da tab_go a tab_DEF
#rende la tabella + leggibile ma con minore informazioni

#tabin="/home/marco/Microarray/csbl/nuovo/a3Ud1/a3Ud1_tab_go"
#tabout="/home/marco/Microarray/csbl/nuovo/a3Ud1/a3Ud1_tab_DEF"

tabgo2tabDEF=function(tabin,tabout) {

tabin=read.table(tabin,sep="\t",stringsAsFactors=F)

righe=gsub("_[A-Z]*","_",rownames(tabin))
righe=unique(righe)
tab=matrix(data=rep("/",5*length(righe)),ncol=5,nrow=length(righe))
colnames(tab)=c("cl","BP","CC","MF","%")

for(i in 1:length(righe)) {
  tab[i,"cl"]=righe[i]
  n=rownames(tabin)[grep(righe[i],rownames(tabin))]
  cat(n,"\n")

  for(j in 1:length(n)) {
    if(length(grep("BP",n[j]))!=0)
      tab[i,"BP"]=tabin[n[j],"desc"]
    if(length(grep("CC",n[j]))!=0)
      tab[i,"CC"]=tabin[n[j],"desc"]
    if(length(grep("MF",n[j]))!=0)
      tab[i,"MF"]=tabin[n[j],"desc"]    
    if(length(grep("TOT",n[j]))!=0 | n[j]=="NN")
      tab[i,"%"]=tabin[n[j],"freq"]

  }

}

a=tab[,"%"]
mode(a)="numeric"
a=round(a/sum(a),2)
tab[,"%"]=a
tab[,"%"]=gsub("0\\.","",tab[,"%"])
tab[,"%"]=gsub("^0","",tab[,"%"])
tab[,"cl"]=gsub("_","",tab[,"cl"])
write.table(tab,file=tabout,quote=FALSE,sep="\t",row.names=FALSE)

}

######################################################################################################



set.prob.table(type="enrichment",filename=file_set)
cat("FATTO set.prob.table!!!!\n deve essere una nuova sessione di R\n")


id_clu=readLines(paste(percorso,nome_ris,"_idclu",sep=""))
id_pul=readLines(paste(percorso,"GO_pul_",nome_ris,sep=""))
tab_go=matrix(data=NA,nrow=length(id_clu)*4,ncol=4)
colnames(tab_go)=c("GO","desc","freq","q.value")

rnames=NULL
for(i in 1:length(id_clu)) {
  rnames=c(rnames,paste("cl",i,"_BP",sep=""))
  rnames=c(rnames,paste("cl",i,"_CC",sep=""))
  rnames=c(rnames,paste("cl",i,"_MF",sep=""))
  rnames=c(rnames,paste("cl",i,"_TOT",sep=""))
}

rownames(tab_go)=rnames

tab_bar=matrix(data=rep(0,4*(length(id_clu)+1)),nrow=4,ncol=length(id_clu)+1)
rownames(tab_bar)=c("BP","CC","MF","tot")
colnames(tab_bar)=c(paste("cl",1:length(id_clu),sep=""),"NN")

system(paste("mkdir ",percorso,"enr",sep=""))


for(i in 1:length(id_clu)) {


  id_sel=unlist(strsplit(id_clu[i]," "))
  ff=file(paste(percorso,"go_temp",sep=""),"w")
  for(j in 1:length(id_sel)) {
    cat(id_pul[grep(id_sel[j],id_pul)],file=ff)
    cat("\n",file=ff)

  }
  close(ff)

  ent=entities.from.text(paste(percorso,"go_temp",sep=""))


  freq=frequencies(ent,p.values=T,fdr=T,normalize.one = FALSE)

  #cat(i,"ok\n")


  tab_bar[4,i]=length(id_sel)


  tab_go[paste("cl",i,"_TOT",sep=""),]=c("TOT","/",length(id_sel),"/")



  if(dim(freq$BP)[1]!=0) {
    tabp=freq$BP
    tabp[,1]=ids.to.descs(tabp[,1])
    tabp=tabp[which(tabp[,6]<0.05),] #qvalue

    if(dim(tabp)[1]!=0) {
      tabw=cbind(tabp,id=rep(NA,dim(tabp)[1]))      
      for(z in 1:dim(tabw)[1]) {
	idr=ent@attrs$key[grep(rownames(tabw)[z],ent@attrs$BP.exp)]
	idr=paste(idr,collapse=", ")
	tabw[z,"id"]=idr
      }
      write.table(tabw,quote=F,file=paste(percorso,"enr/",nome_ris,"_",i,"_freq_GO_BP",sep=""),sep="\t")
      tabp=tabp[order(tabp[,"q.value"],decreasing=FALSE),] #in base a q.value
      tabp[,6]=signif(tabp[,6],2)
      tab_bar[1,i]=tabp[1,2]
      tab_go[paste("cl",i,"_BP",sep=""),]=c(rownames(tabp)[1],tabp[1,1],tabp[1,2],tabp[1,6])
    }
  }

  if(dim(freq$CC)[1]!=0) {
    tabp=freq$CC
    tabp[,1]=ids.to.descs(tabp[,1])
    tabp=tabp[which(tabp[,6]<0.05),] #qvalue

    if(dim(tabp)[1]!=0) {
      tabw=cbind(tabp,id=rep(NA,dim(tabp)[1]))      
      for(z in 1:dim(tabw)[1]) {
	idr=ent@attrs$key[grep(rownames(tabw)[z],ent@attrs$CC.exp)]
	idr=paste(idr,collapse=", ")
	tabw[z,"id"]=idr
      }
      write.table(tabw,quote=F,file=paste(percorso,"enr/",nome_ris,"_",i,"_freq_GO_CC",sep=""),sep="\t")
      tabp=tabp[order(tabp[,"q.value"],decreasing=FALSE),] #in base a q.value
      tabp[,6]=signif(tabp[,6],2)
      tab_bar[2,i]=tabp[1,2]
      tab_go[paste("cl",i,"_CC",sep=""),]=c(rownames(tabp)[1],tabp[1,1],tabp[1,2],tabp[1,6])      
    }
  }

  if(dim(freq$MF)[1]!=0) {
    tabp=freq$MF
    tabp[,1]=ids.to.descs(tabp[,1])
    tabp=tabp[which(tabp[,6]<0.05),] #qvalue
    if(dim(tabp)[1]!=0) {      
      tabw=cbind(tabp,id=rep(NA,dim(tabp)[1]))      
      for(z in 1:dim(tabw)[1]) {
	idr=ent@attrs$key[grep(rownames(tabw)[z],ent@attrs$MF.exp)]
	idr=paste(idr,collapse=", ")
	tabw[z,"id"]=idr
      }
      write.table(tabw,quote=F,file=paste(percorso,"enr/",nome_ris,"_",i,"_freq_GO_MF",sep=""),sep="\t")
      tabp=tabp[order(tabp[,"q.value"],decreasing=FALSE),] #in base a q.value
      tabp[,6]=signif(tabp[,6],2)
      tab_bar[3,i]=tabp[1,2]
      tab_go[paste("cl",i,"_MF",sep=""),]=c(rownames(tabp)[1],tabp[1,1],tabp[1,2],tabp[1,6])
    }
  }
  system(paste("rm ",percorso,"go_temp",sep=""))
}



rem=NULL
for(i in 1:(dim(tab_bar)[2]-1)) {
  if(sum(tab_bar[1:3,i])==0) {
    rem=c(rem,i)
  }
  if(tab_bar[4,i]<n_lim) {
    rem=c(rem,i)
  }
}



rem=unique(rem)

if(length(rem)==(dim(tab_bar)[2]-1)) {
  cat("Rimossi tutti i cluster.\n")
  ff=file(paste(percorso,"cl_anno_",nome_ris,sep=""),"w")
  cat("Sono stati rimossi tutti i clusters.\n",file=ff)
  close(ff)
  return()
}


if(!is.null(rem)) {
  for(i in 1:length(rem)) {
    tab_go=tab_go[-which(rownames(tab_go)==paste("cl",rem[i],"_BP",sep="")),]
    tab_go=tab_go[-which(rownames(tab_go)==paste("cl",rem[i],"_CC",sep="")),]
    tab_go=tab_go[-which(rownames(tab_go)==paste("cl",rem[i],"_MF",sep="")),]
    tab_go=tab_go[-which(rownames(tab_go)==paste("cl",rem[i],"_TOT",sep="")),]
  }
}



if(!is.null(rem))
  tab_bar=tab_bar[,-rem]



rem=NULL
for(i in 1:dim(tab_go)[1]) {
  cat(i,"\n")
  if(is.na(tab_go[i,1]))
    rem=c(rem,i)

}



if(length(rem)==(dim(tab_go)[1]-1)) {
  #da fermarmi o andare avanti
  cat("Rimossi tutti i cluster.\n")
  ff=file(paste(percorso,"cl_anno_",nome_ris,sep=""),"w")
  cat("Sono stati rimossi tutti i clusters.\n",file=ff)
  close(ff)
  return()
}



if(!is.null(rem))
  tab_go=tab_go[-rem,]


n_mac=length(id_pul)-sum(tab_bar[4,])
tab_bar[4,"NN"]=n_mac


tab_go=rbind(tab_go,NN=c("/","/",n_mac,"/"))






tab_go=data.frame(tab_go,stringsAsFactors=FALSE)


write.table(tab_go,quote=FALSE,sep="\t",file=paste(percorso,nome_ris,"_tab_go",sep=""))

tabgo2tabDEF(paste(percorso,nome_ris,"_tab_go",sep=""),paste(percorso,nome_ris,"_tab_DEF",sep=""))

x11(width=14,height=9)
par(mfrow=c(1,2))


tab_bar=cbind(NN=tab_bar[,dim(tab_bar)[2]],tab_bar[,1:(dim(tab_bar)[2]-1)])

barplot(tab_bar, main=nome_ris, xlab="Clusters", ylab="n spots", col=c("lightblue","red","forestgreen","gray"), legend = rownames(tab_bar), beside=TRUE)
textplot(tab_go)
savePlot(filename=paste(percorso,nome_ris,"_bartab.jpg",sep=""),type="jpeg")
graphics.off()

x11(width=13,height=8)
barplot(tab_bar, main=nome_ris, xlab="Clusters", ylab="n spots", col=c("lightblue","red","forestgreen","gray"), legend = rownames(tab_bar), beside=TRUE,cex.axis=2,cex.names=2,cex.main=2)
savePlot(filename=paste(percorso,nome_ris,"_bar.png",sep=""),type="png")
graphics.off()



}
