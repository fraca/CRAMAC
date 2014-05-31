wil_t <-
function(dati,gruppi) {

flag_norm=TRUE

gruppi=as.character(gruppi)
gru=unique(gruppi)
a=dati[which(gruppi==gru[1])]
b=dati[which(gruppi==gru[2])]
sha=shapiro.test(a)
if(sha[["p.value"]]<0.05)
    flag_norm=FALSE
sha=shapiro.test(b)
if(sha[["p.value"]]<0.05)
    flag_norm=FALSE
gruppi=as.factor(gruppi)

flag_omo=TRUE
Ft=var.test(a,b)
if(Ft[["p.value"]]<0.05)
    flag_omo=FALSE

if(flag_norm) {
  if(flag_omo) {
     ris=t.test(a,b, var.equal=TRUE, paired=FALSE)
     print(ris)
  } else {
     ris=t.test(a,b, var.equal=FALSE, paired=FALSE)
     print(ris)
  }

} else {
  ris=suppressWarnings(wilcox.test(a,b))
  print(ris)
}

return(ris$p.value)

}
