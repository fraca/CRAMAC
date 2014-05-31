coo2geni <-
function(coo,db,fout="/tmp/out_coo") {

coo_tot=paste(db,"_coo",sep="")

coo_tot=read.table(coo_tot,sep=";")

coo_tot[,1]=gsub(" seq=cdna","",coo_tot[,1])
rownames(coo_tot)=gsub(">","",coo_tot[,1])
coo_tot[,2]=gsub(" coord=","",coo_tot[,2])
coo_tot[,2]=gsub(":-?1$","",coo_tot[,2])
coo_tot[,1]=gsub(":.+$","",coo_tot[,2])
coo_tot[,3]=gsub("^.+:","",coo_tot[,2])
coo_tot[,2]=gsub("\\.\\..+$","",coo_tot[,3])
coo_tot[,3]=gsub("^.+\\.\\.","",coo_tot[,3])

mode(coo_tot[,2])="numeric"
mode(coo_tot[,3])="numeric"

coo=gsub(",","",coo)
chr=strsplit(coo,":")[[1]][1]
start=strsplit(strsplit(coo,":")[[1]][2],"-")[[1]][1]
end=strsplit(strsplit(coo,":")[[1]][2],"-")[[1]][2]
mode(start)="numeric"
mode(end)="numeric"

ris=coo_tot[which(coo_tot[,1]==chr),]

ris1=ris[which(ris[,2]>start),]
ris1=ris1[which(ris1[,2]<end),]
ris2=ris[which(ris[,2]>start),]
ris2=ris2[which(ris2[,2]<end),]

ris3=rbind(ris1,ris2)
ris3=unique(ris3)

cat(rownames(ris3),sep="\n",file="/tmp/idcdm")
system(paste("fastacmd -d ",db," -i /tmp/idcdm -o ",fout,sep=""))
return(ris3)
}
