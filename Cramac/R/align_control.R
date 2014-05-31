align_control <-
function(primer_L,primer_R,bp_max=10000, e_value=10, lista="/home/marco/maizeseq/ZmB73_RefGen_v2/lista_chr", dir_chr="/home/marco/maizeseq/ZmB73_RefGen_v2/") {
require(seqinr)
require(Biostrings)
primer_control=function(primer, evalue, lista, dir_chr) {

write.fasta(primer,names="bla",file.out="/tmp/seq.fasta")
lista=readLines(lista)

ff=file("/tmp/blast_tab","w")
for(i in 1:length(lista)) {
    testo_blast=paste("blast2 -i /tmp/seq.fasta -d ",dir_chr,lista[i],"/",lista[i], " -e ",evalue," -p blastn -m  8 -o /tmp/tab_ris",sep="")
    system(testo_blast)
    a=readLines("/tmp/tab_ris")
    #cat(i,"\n")
    if(length(a)!=0) {
      for(j in 1:length(a)) {
	cat(lista[i],"\t",a[j],"\n",sep="",file=ff)

      }
    }
}
close(ff)
tab=read.table("/tmp/blast_tab",sep="\t",stringsAsFactors=FALSE)

for(i in 1:dim(tab)[1]) {
  if(tab[i,10]<tab[i,11]) {
    tab[i,2]="L"
    tab[i,3]=tab[i,11]
  } else {
    tab[i,2]="R"
    tab[i,3]=tab[i,10]
  }

}
tab=tab[,1:3]

return(tab)
}

s1=DNAString(primer_L)
s2=DNAString(primer_R)
mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
ali=pairwiseAlignment(s1, s2, substitutionMatrix = mat,  gapOpening = -1, gapExtension = -1)
print(ali)
tabL=primer_control(primer_L, evalue=e_value, lista=lista, dir_chr=dir_chr)
tabR=primer_control(primer_R,evalue=e_value,lista=lista,dir_chr=dir_chr)
cat("\n","primerL n match:",dim(tabL)[1],"  primerR match:",dim(tabR)[1],"\n")
chr_sel=intersect(unique(tabL[,1]),unique(tabR[,1]))

for(i in 1:length(chr_sel)) {
  selL=tabL[which(tabL[,1]==chr_sel[i]),]
  selR=tabR[which(tabR[,1]==chr_sel[i]),]
  mode(selL[,3])="numeric"
  mode(selR[,3])="numeric"

  for(j in 1:dim(selL)[1]) {
    for(z in 1:dim(selR)[1]) {
      if(selL[j,2]!=selR[z,2]) {
	bp=abs(selL[j,3]-selR[z,3])
	if(bp<10000)
	  cat(selL[j,1],selL[j,2],selL[j,3],selR[z,2],selR[z,3],bp,"\n")
      }
    }
  }

}

}
