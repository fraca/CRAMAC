primer_s <-
function(seq,temp=60,win_min=19,win_max=21,n_primer=5,fin=TRUE) {
#  win_min=19
#  win_max=23
#  n_primer=5
# temp=60
# # seq="TTGTTTTAACATAAAACAATTTGTATTCCTAATAATAGAATATTAAAAAAGAAAAGCTACCAATTAGCTCTTTCATGATTACCTCATCAAATAAGTTTTAATGTTACCATCTTAATTTATATCAGGCT"
# seq="AGGGATATTCTTTGTCTGGACCAGAAAGATTGGATTGTACTCGATCGGGACGCTGGACAGACTCCCCACCAATGTGTGAAG"
#  fin=TRUE
#  n_primer=5

seq=toupper(seq)
if(gsub("[A,T,C,G]","",seq)!="")
  stop("Error type unknown.\n")
seq=unlist(strsplit(seq,""))

win_tot=seq(win_min,win_max,1) ##scegli bp primer
n_temp=NULL
primer=NULL
for(j in 1:length(win_tot)) {
  win=win_tot[j]
  inizi=seq(1,(length(seq)-win),1)
  fini=seq(win,(length(seq)-1),1)
  #length(inizi)
  #length(fini)

  for(i in 1:length(inizi)) {
    #cat(i,"\n")
    flag=TRUE
    if(fin) {
      if(seq[inizi[i]]=="A" | seq[inizi[i]]=="T" | seq[fini[i]]=="A" | seq[fini[i]]=="T")
	flag=FALSE
    }

    oligo=0
    n1=seq[inizi[i]]
    n2=seq[(inizi[i]+1)]
    z=inizi[i]+2


    while (z<=fini[i] & flag==TRUE) {
      #cat(z,fini[i],"\n")
      if( (z+2)<=fini[i] & (z-3)>=inizi[i]) { #palindromi di 6
	#cat("entro\n")
	n3=seq[z]
	n4=seq[z+1]
	n5=seq[z+2]
	n_3=seq[z-1]
	n_4=seq[z-2]
	n_5=seq[z-3]
	#cat(paste(n3,n4,n5,sep=""),toupper(paste(revcomp(n_5),revcomp(n_4),revcomp(n_3),sep="")),"\n")
	if(paste(n3,n4,n5,sep="")==toupper(revcomp(paste(n_5,n_4,n_3,sep=""))))
	  flag=FALSE
      }
      z=z+1
    }
    pri_a=seq[inizi[i]:fini[i]]
    pri=paste(pri_a,collapse="")
    if(length(grep("AAAA",pri))!=0)
      flag=FALSE
    if(length(grep("TTTT",pri))!=0)
      flag=FALSE
    if(length(grep("CCCC",pri))!=0)
      flag=FALSE
    if(length(grep("GGGG",pri))!=0)
      flag=FALSE


    if(flag) {   
    #n_GC=which(pri_a=="C" | pri_a=="G")
    #n_GC=length(n_GC)
    #n_GC=c(n_GC,nchar(gsub("[A,T]","",pri)))
    #n_temp=c(n_temp,(64.9 + 41*(nchar(gsub("[A,T]","",pri))-16.4)/nchar(pri)))
    n_temp=c(n_temp,(69.3 +( 41*nchar(gsub("[A,T]","",pri))/ nchar(pri)-(650 / nchar(pri))))) #presa da eurofins
    #cat(n_temp,"\n")
    #per_GC=c(per_GC,(length(n_GC)/win))
    primer=c(primer,pri)
    }
  }

}

o_GC=abs(n_temp-temp)
o_GC=order(o_GC)



if(length(o_GC)< n_primer) 
  n_primer=length(o_GC)

cat("\n\nThe primers:\n\n")
for(i in 1:n_primer) { ##stampa i migliori

  cat(primer[o_GC[i]],"\n")
  cat(n_temp[o_GC[i]],"\n")
  cat(nchar(primer[o_GC[i]]),"\n\n")

}

}
