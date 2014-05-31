t_mel <-
function(primer) {
primer=toupper(primer)
if(gsub("[A,T,C,G]","",primer)!="")
  stop("Error type unknown.\n")
#a=nchar(gsub("[A,T]","",primer))/nchar(primer)
#a=(64.9 + 41*(nchar(gsub("[A,T]","",primer))-16.4)/nchar(primer)) presa da un sito
a= 69.3 +( 41*nchar(gsub("[A,T]","",primer))/ nchar(primer)-(650 / nchar(primer)) ) #presa da eurofins
cat(a,"\n")
return(a)
}
