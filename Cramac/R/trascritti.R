trascritti <-
function(file_gen,file_html) {

gb2fasta(file_gen,"/tmp/ensembl.fasta")

seq=read.fasta("/tmp/ensembl.fasta") 
# seq=readLines("/tmp/ensembl.fasta") ## se read.fasta non funziona  (errore strano s2c) prova questo probabilmente e' da rinstallare R
# seq=seq[-1]
# seq=list(unlist(strsplit(seq,"")))

gen=readLines(file_gen)

pos_mRNA=grep("mRNA",gen)
cat("\nCi sono",length(pos_mRNA)," trascritti\n")

tra=NULL
nomi_tra=NULL
for(i in 1:length(pos_mRNA)) {
  flag=TRUE
  flag2=TRUE
  bps=NULL
  n_riga=pos_mRNA[i]
  while(flag) {    
    if(length(grep("/note=",gen[n_riga]))==0) {
      #if(length(grep("\\)",gen[n_riga]))!=0)
      if(length(grep("/",gen[n_riga]))!=0)
	flag2=FALSE
      if(flag2) {
	riga=gsub("join\\(","",gen[n_riga])
	riga=gsub(",$","",riga)
	riga=gsub("[\\) ]","",riga)
	riga=gsub("mRNA","",riga)
	#cat(riga,"\n")
	coppie=unlist(strsplit(riga,","))	
	for(z in 1:length(coppie)) {
	  cop=unlist(strsplit(coppie[z],"\\.\\."))
	  mode(cop)="numeric"
	  #cat(pos_mRNA[i],cop[1],cop[2],"\n")
	  bps=c(bps,cop[1]:cop[2])
	}
      }
      
    } else {
      nome=gsub("/note=\"transcript_id=","",gen[n_riga])
      nome=gsub("/note=\"identifier=","",nome)
      nome=gsub("[\", ]","",nome)      
      nomi_tra=c(nomi_tra,nome)
      flag=FALSE
    }
    n_riga=n_riga+1

  }
  tra[[i]]=bps
  #cat("\n\n")
}
#cat(length(tra),length(nomi_tra),"\n")
names(tra)=nomi_tra

mat=matrix(data=0,nrow=length(seq[[1]]),ncol=length(tra))
righe=NULL
for(i in 1:dim(mat)[1]) {
  for(j in 1:dim(mat)[2]) {
    if(length(which(tra[[j]]==i))!=0)
      mat[i,j]=1
  }
  righe=c(righe,paste(mat[i,],collapse=""))

}


un_righe=unique(righe)
col0=paste(rep(0,length(tra)),collapse="")
un_righe=un_righe[-which(un_righe==col0)]

colori=c("#0000FF","#DC143C","#008800","#FFA500","#B8860B","#8A2BE2","#A52A2A","#5F9EA0","#7FFF00","#D2691E","#FF7F50","#00008B", "#FF1493","#800000","#FFFF00","#9ACD32","#FF4500","#98FB98","#A9A9A9","#DEB887")
if(length(un_righe)>length(colori)) 
  stop("ci sono piu\' di 20 combinazioni\naumentare colori.....\n")
colori=colori[1:length(un_righe)]
names(colori)=un_righe

ff=file(file_html,"w")

cat("<HTML>\n\n",file=ff)
cat("<PRE><TT>\n",file=ff)
cat("<B>legenda</B>\n",file=ff)
cat("introne in tutti i trascritti\n",file=ff)

for(i in 1:length(colori)) {
  ri=unlist(strsplit(un_righe[i],""))
  ri=gsub("1",TRUE,ri)
  ri=gsub("0",FALSE,ri)
  mode(ri)="logical"
  cat("<B><span style='color:",colori[i],";'>",i," ",paste(nomi_tra[ri],collapse=" "),"</span></B>\n",sep="",file=ff)
}
cat("\n\n<B>>sequenza</B>\n",file=ff)

riga_p=righe[1]
if(riga_p!=col0)
  cat("<B><span style='color:",colori[riga_p],";'>",sep="",file=ff)

for(i in 1:length(righe)) {

  c1=i/60
  c2=round(c1)
  if(c1==c2)
    cat("\n",file=ff)
  if(righe[i]!=riga_p) {
    cat("</span></B>",file=ff)
    riga_p=righe[i]
    if(righe[i]!=col0)
      cat("<B><span style='color:",colori[righe[i]],";'>",sep="",file=ff)
  }
  cat(toupper(seq[[1]][i]),file=ff)
}

if(righe[i]!=col0)
  cat("</span></B>",file=ff)

close(ff)

}
