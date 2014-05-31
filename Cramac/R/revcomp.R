revcomp <-
function(seq,lowercase=FALSE) {
  seq=toupper(seq)
  seq=unlist(strsplit(seq,""))
  seq2=NULL
  for(i in length(seq):1) {
    flag=TRUE

    if(seq[i]=="A" | seq[i]=="a") {
      seq2=c(seq2,"T")
      flag=FALSE
    }
    if(seq[i]=="T"| seq[i]=="t") {
      seq2=c(seq2,"A")
      flag=FALSE
    }

    if(seq[i]=="G"| seq[i]=="g") {
      seq2=c(seq2,"C")
      flag=FALSE
    }
    if(seq[i]=="C"| seq[i]=="c") {
      seq2=c(seq2,"G")
      flag=FALSE
    }


    if(flag)
      stop("Error type unknown.",seq[i],"\n")
  }
  seq2=paste(seq2,collapse="")
  if(lowercase)
    seq2=tolower(seq2)
  return(seq2)
}
