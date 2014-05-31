csbl1 <-
function(id,file_set,nome_GO,nome_ris,percorso,cutoff=0.5,metric="Lin",onto=NULL) {




set.prob.table(type="similarity",filename=file_set)
cat("FATTO set.prob.table!!!!\n deve essere una nuova sessione di R\n")

#leggo GO e tolgo MZ che non hanno GO

a=readLines(nome_GO)
a_sp=strsplit(a," ")
ff=file(paste(percorso,"GO_pul_",nome_ris,sep=""),"w")

#usa solo un ontologia
rem_a=NULL
if(!is.null(onto)) {
  for(i in 1:length(a_sp)) {
    rem_j=NULL
    for(j in 2:length(a_sp[[i]])) {
      if(ids.to.ontologies(a_sp[[i]][j])!=onto) 
	rem_j=c(rem_j,j)
    }
    if(!is.null(rem_j))
      a_sp[[i]]=a_sp[[i]][-rem_j]
    if(length(a_sp[[i]])==1) {
      rem_a=c(rem_a,i)
    } else {
      a[i]=paste(a_sp[[i]],collapse=" ")
    }
  }
}

if(!is.null(rem_a))
  a=a[-rem_a]
a_sp=strsplit(a," ")

for(i in 1:length(a)) {
  wh=which(id==a_sp[[i]][1])
  if(length(wh)!=0) {
    if(length(a_sp[[i]])!=1)
      cat(a[i],"\n",sep="",file=ff)
  }
}
close(ff)

#nome geni tab1 e tab2
tab_h=matrix(data=rnorm(length(id)*2),nrow=length(id),ncol=2)
rownames(tab_h)=id
heatb=go.heatmap(tab_h,paste(percorso,"GO_pul_",nome_ris,sep=""),grayscale=FALSE,go.cut=cutoff,use.inform=FALSE,metric=metric)
graphics.off()

ff=file(paste(percorso,nome_ris,"_idclu",sep=""),"w")
for(j in 1:length(heatb$members)) {
  cat(heatb$members[[j]],sep=" ",file=ff)
  cat("\n",file=ff)
}

close(ff)
}
