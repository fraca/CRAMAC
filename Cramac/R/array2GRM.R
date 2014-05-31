array2GRM <-
function(id_array,tab_conf,tab_array=NULL,transcript=FALSE,score_length=NULL,file_gramene="/tmp/X_gramene") {

tab_c=read.table(tab_conf,header=TRUE,stringsAsFactors=FALSE)

if(!is.null(score_length))
  tab_c=tab_c[which(tab_c[,3]>score_length),]

tab_ris=tab_c[1,]
rimossi=NULL
for(i in 1:length(id_array)) {
  nn=which(tab_c[,1]==id_array[i])
  if(length(nn)!=0) {
    for(j in 1:length(nn)) {
      tab_ris=rbind(tab_ris,tab_c[nn[j],])
    }
  } else {
    rimossi=c(rimossi,id_array[i])
  }
}
cat("Sono stati rimossi",length(rimossi),"su",length(id_array),"id array che non sono presenti in tabella (/tmp/rimossi)\n")
cat(rimossi,sep="\n",file="/tmp/rimossi")
tab_ris=tab_ris[-1,]

if(!transcript) {
  tab_ris[,2]=gsub("_[A-Z,0-9]+","",tab_ris[,2])
  GRM_u=strsplit(unique(paste(tab_ris[[1]],tab_ris[[2]],sep="-")),"-")

  tab_ris2=tab_ris[1,]
  for(i in 1:length(GRM_u)) {
    righe=intersect(which(tab_ris[,1]==GRM_u[[i]][1]),which(tab_ris[,2]==GRM_u[[i]][2]))
    tab_ris2=rbind(tab_ris2,tab_ris[righe[which(tab_ris[righe,"score_length"]==max(tab_ris[righe,"score_length"]))[1]],])
  

  }
  tab_ris2=tab_ris2[-1,]
  tab_id=tab_ris2
} else {
  tab_id=tab_ris
}



if(!is.null(tab_array)) {
  flag_v=FALSE
  if(is.vector(tab_array)) {
    names(tab_array)=id_array
    tab2=cbind(tab_id,val=tab_array[tab_id[,1]])
    flag_v=TRUE
  } else {
    tab2=cbind(tab_id,tab_array[tab_id[,1],])
  }
  GRM_u=unique(tab2[,2])
  tab3=tab2[1,]
  tab4=tab2[1,]
  #id_fatti=NULL
  for(i in 1:length(GRM_u)) {
    nn=which(tab2[,2]==GRM_u[i])
    if(length(nn)==1) {
      tab3=rbind(tab3,tab2[nn,])
    } else {
      ##da sistemare
      ida=tab2[nn,1]
      tab2[nn[1],1]=paste(tab2[nn,1],round(tab2[nn,3],2),sep="_",collapse="-")
      if(flag_v)
	tab2[nn[1],5]=mean(tab2[nn,5])
      else
	tab2[nn[1],c(5:dim(tab2)[2])]=apply(tab2[nn,c(5:dim(tab2)[2])],2,mean)

      tab4=rbind(tab4,tab2[nn[1],])
      tab2[nn,2]=rep("bla",length(nn)) #non metto roba in piu'
      cat(GRM_u[i],"ha piu\' id_array",ida,"(fatta media)\n")
    }
  }
  tab3=tab3[-1,]
  tab4=tab4[-1,]

  tab_fin=rbind(tab3,tab4)
  tab_fin=tab_fin[,c(1,2,5:dim(tab_fin)[2])]
  ff=file(file_gramene,"w")

  write.table(tab3[,c(2,5:dim(tab3)[2])],sep="\t",quote=FALSE,row.names=FALSE,file=ff)
  cat("\n####################################\nquesti hanno piu\' id_array\n##############################\n",file=ff)
  write.table(tab4[,c(2,5:dim(tab4)[2])],sep="\t",quote=FALSE,row.names=FALSE,file=ff)
  close(ff)
  ris_tot=list(tab_id=tab_id,tab_array=tab_fin)

} else {
 ris_tot=tab_id
 
}
return(ris_tot)

}
