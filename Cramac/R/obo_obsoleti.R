obo_obsoleti <-
function(obo_tot,nome) {

obo_tot=readLines(obo_tot)

sel=grep("is_obsolete: true",obo_tot)

ff=file(nome,"w")
cat("GO_obsoleti\tGO_nuovi\n",file=ff)
for(i in 1:length(sel)) {
#for(i in 1:3) {
  flag=TRUE
  j=1
  while(flag) {
    if(length(grep("^id: ",obo_tot[sel[i]-j]))!=0) {
      flag=FALSE
      cat(gsub("^id: ","",obo_tot[sel[i]-j]),"\t",sep="",file=ff)
    }
    j=j+1
    #cat(i,j,"\n")
  }
  flag=TRUE
  flag2=FALSE
  j=1
  while(flag) {
    if(length(grep("^replaced_by: ",obo_tot[sel[i]+j]))!=0) {
      if(flag2) {
	cat(gsub("^replaced_by:","",obo_tot[sel[i]+j]),file=ff)
      } else {
	flag2=TRUE
	cat(gsub("^replaced_by: ","",obo_tot[sel[i]+j]),file=ff)
      }
    }
    if(length(grep("^consider: ",obo_tot[sel[i]+j]))!=0) {
      if(flag2) {
	cat(gsub("^consider:","",obo_tot[sel[i]+j]),file=ff)
      } else {
	flag2=TRUE
	cat(gsub("^consider: ","",obo_tot[sel[i]+j]),file=ff)
      }
    }
    if(length(grep("^\\[Term\\]",obo_tot[sel[i]+j]))!=0) {
      flag=FALSE
      if(flag2) {
	cat("\n",file=ff)
      } else {
	cat("NA\n",file=ff)
      }
    }
    #cat(sel[i],j,"\n")
    j=j+1
  }

}

sel=grep("^alt_id: ",obo_tot)


for(i in 1:length(sel)) {
  cat(gsub("^alt_id: ","",obo_tot[sel[i]]),"\t",sep="",file=ff)
  flag=TRUE
  j=1
  while(flag) {
    if(length(grep("^id: ",obo_tot[sel[i]-j]))!=0) {
      flag=FALSE
      cat(gsub("^id: ","",obo_tot[sel[i]-j]),"\n",sep="",file=ff)
    }
    j=j+1
    cat(i,j,"\n")
  }


}


close(ff)


}
