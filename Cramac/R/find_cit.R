find_cit <-
function(file_in,file_out) {

a=readLines(file_in)

ris=NULL
for(i in 1:length(a)) {

if(substr(a[i],1,1)=="(" & substr(a[i],1,1)=="(")
  stop("gjgjhgjkjl")

  if(a[i]!="") {
    riga=unlist(strsplit(a[i],""))
    ent=""
    for(j in 1:length(riga)) {
      if(riga[j]=="(") {
	s=(j-30)
	e=j
	flag=TRUE
	while(flag) {
	  if(riga[e]==")")
	    flag=FALSE
	  else
	    e=e+1

	}

	if(s<=0)
	  s=1
	ris=c(ris,paste(riga[s:e],sep="",collapse=""))
	} 

    }



  }

  #cat(i,"\n")
}
ff=file(file_out,"w")
cat(ris,sep="\n",file=ff)
close(ff)

}
