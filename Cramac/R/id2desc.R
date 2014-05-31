id2desc <-
function(id,tipo,fout="/tmp/desc_out") {

ff=file(fout,"w")
desc_tot=NULL
if(tipo=="array") {
  anno=readLines("/home/marco/maizeseq/maizearray.v4.annotation.txt")
  for(i in 1:length(id)) {
    desc=anno[grep(id[i],anno)]
    for(j in 1:length(desc)) {
      cat(desc[j],"\n",file=ff)

    }

    desc_tot=c(desc_tot,desc)
  }

}
if(tipo=="affy") {
 anno=readLines("/home/marco/maizeseq/Maize.na32.annot.csv")
  for(i in 1:length(id)) {
    desc=anno[grep(id[i],anno)]
    for(j in 1:length(desc)) {
      cat(desc[j],"\n",file=ff)

    }

    desc_tot=c(desc_tot,desc)
  }

}
if(tipo=="GRM") {
  anno=readLines("/home/marco/maizeseq/display.txt")
  anno2=readLines("/home/marco/maizeseq/xref.txt")

  idp=gsub("_.+$","",id)
  idp=unique(idp)


  for(i in 1:length(idp)) {
    cat(idp[i],"\n\n",file=ff)
  
    cat("trascritti:",id[grep(idp[i],id)],"\n\n",file=ff)

    desc=anno[grep(idp[i],anno)]
    desc=c(desc,anno2[grep(idp[i],anno2)])

    for(j in 1:length(desc)) {
      cat(desc[j],"\n",file=ff)

    }
    cat("\n\n",file=ff)
    desc_tot=c(desc_tot,desc)
  }

}

close(ff)
return(desc_tot)
}
