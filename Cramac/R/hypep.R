hypep <-
function(ids,db) {

Doo_tot=NULL
OHM_tot=NULL
MW_tot=NULL
residue_tot=NULL
inc_tot=NULL
isoel_tot=NULL
id_tot=NULL
pro_res=matrix(data=NA,ncol=9)
colnames(pro_res)=c("Tiny", "Small", "Alip", "Arom", "Nonp", "Pol", "Char", "Bas", "Acidic")

for(i in 1:length(ids)) {
  id=ids[i]
  id=gsub("_FGT","_FGP",id)
  id=gsub("_T","_P",id)
  id_tot=c(id_tot,id)
  system(paste("fastacmd -s ",id," -d ",db," -o /tmp/seq.fasta",sep=""))
  system("pepinfo -seq /tmp/seq.fasta -out /tmp/pep_ris -goutfile /tmp/pep -graph cps -generalplot no")
  system("pepstats -seq /tmp/seq.fasta -out /tmp/pepstat_ris")
  ##per stat varie
  sta=readLines("/tmp/pepstat_ris")
  MW=gsub("[M,a-z,=, ,\t]+","",sta[3])
  residue=unlist(strsplit(MW,"R"))[2]
  MW=unlist(strsplit(MW,"R"))[1]
  mode(MW)="numeric"
  mode(residue)="numeric"
  MW_tot=c(MW_tot,MW)
  residue_tot=c(residue_tot,residue)
  inc=unlist(strsplit(sta[8]," = "))[2]
  mode(inc)="numeric"
  inc_tot=c(inc_tot,inc)
  isoel=unlist(strsplit(sta[5]," = "))[2]
  mode(isoel)="numeric"
  isoel_tot=c(isoel_tot,isoel)    
  ri=NULL
  for(j in 39:47) {
    ca=gsub("\t\t\t","\t",sta[j])
    ca=gsub("\t\t","\t",sta[j])
    ca=unlist(strsplit(ca,"\t"))[4]
    mode(ca)="numeric"
    ri=c(ri,ca)
  }
  pro_res=rbind(pro_res,ri)
  ##per Hydro
  hy=readLines("/tmp/pep_ris")
  hy=hy[-grep("^$",hy)]
  hy=hy[-grep("Position",hy)]
  flag1=FALSE
  flag2=FALSE
  Doo=0
  OHM=0
  for(j in 1:length(hy)) {
    if(length(grep("Consensus",hy[j]))!=0)
	flag2=FALSE
    if(flag2) {
      ri=unlist(strsplit(hy[j],"[ ]+"))[4]
      mode(ri)="numeric"
      OHM=OHM+ri
    }
    if(length(grep("OHM",hy[j]))!=0) {
      flag1=FALSE
      flag2=TRUE
    }
    if(flag1) {
      ri=unlist(strsplit(hy[j],"[ ]+"))[4]
      mode(ri)="numeric"
      Doo=Doo+ri
    }
    if(length(grep("Doolittle",hy[j]))!=0)
      flag1=TRUE
  }
  Doo=Doo/residue
  OHM=OHM/residue
  Doo_tot=c(Doo_tot,Doo)
  OHM_tot=c(OHM_tot,OHM)
}

pro_res=pro_res[-1,]
rownames(pro_res)=id_tot
ris=list(id=id_tot,HY_Doo=Doo_tot,HY_OHM=OHM_tot,MW=MW_tot,n_AA=residue_tot,inc_bodies=inc_tot,tab_AA=pro_res)
return(ris)
}
