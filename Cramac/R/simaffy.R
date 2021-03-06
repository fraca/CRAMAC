simaffy <-
function(dir,conf,fold_change=0.58,p_value=0.05,nome_covdesc="covdesc.txt",heat_flag=FALSE,graph_flag=FALSE) {

require(simpleaffy)

setwd(dir)

x <- read.affy(nome_covdesc) # Reads data in working directory.
x.mas5 <- call.exprs(x,"mas5") # Calculates expression values with MAS 5.0 method which is required for the next step!
qc <- qc(x,x.mas5)

plot(qc)

savePlot(filename="Quality_control.jpg",type="jpeg")
graphics.off()
flag=TRUE
id_sel=NULL
for(i in 1:length(conf)) {

  comp <- pairwise.comparison(x.mas5,"treatment",conf[[i]],x)

  if(flag) {
    tab_tot=cbind(comp@fc,comp@tt)
    colnames(tab_tot)=c(paste(paste(conf[[i]],collapse="_"),"fc",sep="_"), paste(paste(conf[[i]],collapse="_"),"tt",sep="_"))
    flag=FALSE
  } else {
    if(!unique(rownames(tab_tot)==names(comp@fc)))
      stop("ordine array diverso!!!!\n")
    no_old=colnames(tab_tot)
    tab_tot=cbind(tab_tot,comp@fc,comp@tt)
    colnames(tab_tot)=c(no_old,paste(paste(conf[[i]],collapse="_"),"fc",sep="_"), paste(paste(conf[[i]],collapse="_"),"tt",sep="_"))
  }
  if(graph_flag) {
   plot(comp,type="scatter")
   savePlot(filename=paste("plot_",conf[[i]][1],"_",conf[[i]][2],"_scatter.jpg",sep=""),type="jpeg")
   graphics.off()
 
   plot(comp,type="ma")
   savePlot(filename=paste("plot_",conf[[i]][1],"_",conf[[i]][2],"_ma.jpg",sep=""),type="jpeg")
   graphics.off()
 
   plot(comp,type="volcano")
   savePlot(filename=paste("plot_",conf[[i]][1],"_",conf[[i]][2],"_volcano.jpg",sep=""),type="jpeg")
   graphics.off()
  }
  
  fil <- pairwise.filter(comp,min.present.no="all",present.by.group=T,fc=fold_change,tt= p_value)

#   if(heat_flag) { #vecchio non usare
#     hmap.pc(fil,x.mas5)
#     savePlot(filename=paste("hetmap_",conf[[i]][1],"_",conf[[i]][2],".jpg",sep=""),type="jpeg")
#     graphics.off()
#   }

  fc=fil@fc
  tt=fil@tt
  a=cbind(fc,tt)
  a=a[order(a[,1]),]
  id_sel=c(id_sel,rownames(a))
  write.table(a,quote=FALSE,sep="\t",file=paste("ris_",conf[[i]][1],"_",conf[[i]][2],".txt",sep=""))
  rem=NULL


}

#vecchio metodo tiene dentro anche AFF e altre schifezze
# tt_col=grep("_tt$",colnames(tab_tot))
# fc_col=grep("_fc$",colnames(tab_tot))
# #cat(dim(tab_tot)[1],"\n")
# for(i in 1:dim(tab_tot)[1]) {
#   if(length(unique(tab_tot[i,tt_col]<p_value))==1)
#    if(!unique(tab_tot[i,tt_col]<p_value)) #tolgo se maggiore di pvalue
#       rem=c(rem,i)
#   if(length(unique(abs(tab_tot[i,fc_col])<fold_change))==1)
#     if(unique(abs(tab_tot[i,fc_col])<fold_change)) #tolgo se minore di fold change
#       rem=c(rem,i)
# }
# rem=unique(rem)
# tab_tot=tab_tot[-rem,]

tab_tot=tab_tot[unique(id_sel),]
write.table(tab_tot,quote=FALSE,sep="\t",file="tab_tot")

if(heat_flag) {
  cat(dim(tab_tot)[1],"sopra mille e\' molto lento...\n")
  fc_col=grep("_fc$",colnames(tab_tot))
  tab_fc=tab_tot[,fc_col]
  colnames(tab_fc)=gsub("_fc$","",colnames(tab_fc))
  colnames(tab_fc)=gsub("_","",colnames(tab_fc))
  require(gplots)
  x11(width=14,height=9)
  hea=heatmap.2(tab_fc, col=redgreen(75), scale="row", key=TRUE, symkey=FALSE, density.info="none", trace="none",cexCol=1, cexRow=0.5,margins=c(3,5),ylab="SPOTS")
  save(tab_fc,hea,file=paste(dir,"hea.Rdata"))
  savePlot(filename=paste(dir,"heatmap_tot.png",sep=""),type="png")
  graphics.off()
}
#cat(dim(tab_tot)[1],"\n")
return(tab_tot)

}
