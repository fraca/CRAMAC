anva <-
function(dati,gruppi,p_value=0.05) {

require(agricolae)

flag_norm=TRUE
gruppi=as.character(gruppi)
gru=unique(gruppi)
if(length(gru)<=2) 
  stop("Utilizza wil_t.R solo 2 livelli della varabile indipendente\n")
ris=matrix(data=NA,nrow=length(gru),ncol=3)
rownames(ris)=unique(gruppi)
colnames(ris)=c("mean","sd","Comparison")
for(i in 1:length(gru)) {
  dat=dati[which(gruppi==gru[i])]
  ris[gru[i],"mean"]=round(mean(dat),2)
  ris[gru[i],"sd"]=round(sd(dat),2)
  sha=shapiro.test(dat)
  if(sha[["p.value"]]<0.05) {
    flag_norm=FALSE
  }
}

ris=ris[order(ris[,1]),]

gruppi=as.factor(gruppi)


flag_omo=TRUE


bar=bartlett.test(dati, gruppi)
    
if(bar[["p.value"]]<0.05) {
  flag_omo=FALSE
}

fli=fligner.test(dati, gruppi)
if(fli[["p.value"]]<0.05) {
  flag_omo=FALSE
}

if(flag_norm) {
  if(flag_omo) {
    modello = lm(formula = dati ~ gruppi)
    ano=anova (modello)
    print(ano)

    if(ano[["Pr(>F)"]][1]<p_value) {
      df<-df.residual(modello)
      MSerror<-deviance(modello)/df
      co <- LSD.test(dati, gruppi, df, MSerror, group=TRUE)
      dun=gsub(" ","",as.character(co[["M"]]))
      names(dun)=gsub(" ","",as.character(co[["trt"]]))
      ris[,"Comparison"]=dun[rownames(ris)]
    }
  } else {
    data=cbind(gruppi,dati)
    wel=oneway.test(dati ~ gruppi, data=data, var.equal=FALSE)
    print(wel)
    if(wel[["p.value"]]<p_value) {
      model <- lm(dati ~ gruppi)  
      df<-df.residual(model)  
      MSerror<-deviance(model)/df  
      Fc<-anova(model)[1,4]
      co <- waller.test(dati, gruppi, df, MSerror, Fc, group=TRUE)
      dun=gsub(" ","",as.character(co[["M"]]))
      names(dun)=gsub(" ","",as.character(co[["trt"]]))
      ris[,"Comparison"]=dun[rownames(ris)]
    }
  }
} else {
  if(flag_omo) {
    co=kruskal(dati,gruppi)
    dun=gsub(" ","",as.character(co[["M"]]))
    if(length(unique(dun))!=1) {
      names(dun)=gsub(" ","",as.character(co[["trt"]]))
      ris[,"Comparison"]=dun[rownames(ris)]
    }
  } else {
    nn=length(which(gruppi==levels(gruppi)[1]))
    judge=rep(1:nn,length(levels(gruppi)))
    co<-friedman(judge,gruppi,dati)
    dun=gsub(" ","",as.character(co[["M"]]))
    if(length(unique(dun))!=1) {
      names(dun)=gsub(" ","",as.character(co[["trt"]]))
      ris[,"Comparison"]=dun[rownames(ris)]
    }
  }
}
return(ris)
}
