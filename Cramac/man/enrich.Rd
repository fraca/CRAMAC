\name{enrich}
\alias{enrich}

\title{
enrich
}
\description{
The function does enrichment analysis through csbl.go packages.
}
\usage{
enrich(id, nome_ris, percorso, file_set, nome_GO)
}

\arguments{
  \item{id}{
Vector with the selected ids.
}
  \item{nome_ris}{
Filename of output.
}
  \item{percorso}{
Directory where to put the output.
}
  \item{file_set}{
Filename of enrichment probability table.
}
  \item{nome_GO}{
File with ids and GO, one id per line.
}
}
\details{
GENE ONTOLOGY ANALYSIS: \cr
1. obo_obsoleti find obsolete GO term in obo tot file. \cr
2. rem_obo remove obsolete GO term for GO dataset. \cr
3. redGO reduce GO dataset by GOslim. \cr
4. prob_table calculate similarity and enrichment probability tables. \cr
5. enrich  enrichment analysis. \cr
6. csbl1 clustering based on the GO similarity matrix \cr
7. csbl2 clustering based on the GO similarity matrix \cr
}

\value{
For each ontology: \cr \cr
Table_freq (goid,desc,freq,proportion,p.value,priori,q.value)\cr
table with all GO terms. \cr
goid GO term. \cr
desc description of GO term. \cr
freq number of ids with GO term. \cr
proportion proportion of ids with GO term. \cr
p.value p value using Fisher's Exact Test (from csbl.go) \cr
priori value from enrichment probability table. \cr
q.value false discovery rate (from csbl.go) \cr \cr
Table_sel (goid,query,ref,desc,p_value,q_value,ids) \cr 
table with selected GO terms (query>ref, p_value < 0.05, q_value < 0.05). \cr
goid GO term. \cr
query proportion of ids with GO term. \cr
ref proportion of ids with GO term in reference. \cr
desc description of GO term. \cr
p.value p value using Fisher's Exact Test (from csbl.go) \cr
q.value false discovery rate (from csbl.go) \cr
ids ids with the GO term. \cr \cr
Bargraph with selected GO terms. 
}
\author{
Marco Fracassetti
}

\seealso{
\code{\link{obo_obsoleti}} \cr
\code{\link{rem_obo}} \cr
\code{\link{redGO}} \cr
\code{\link{prob_table}} \cr
\code{\link{enrich}} \cr
\code{\link{csbl1}} \cr
\code{\link{csbl2}} \cr
}

\examples{
# id=readLines("/home/marco/biosegen_array/definitivo/setdiff_GRM")
# nome_ris="setdiff"
# percorso="/home/marco/biosegen_array/definitivo/csbl_pro/"
# file_set="/home/marco/funzioni_R/agriGO/prob_tab_sim_slimGRM_def"
# nome_GO="/home/marco/funzioni_R/agriGO/slimGO_GRM_csbl_pul"
# enrich(id,nome_ris,percorso,file_set,nome_GO)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
