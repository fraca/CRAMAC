\name{csbl2}
\alias{csbl2}

\title{
csbl2
}
\description{
 Agglomerative hierarchical clustering based on the GO similarity matrix. \cr
 For each cluster a GO term enrichment analysis was performed, \cr
 for each of the three ontology class (BP, CC, MF), \cr
 the GO term with the maximum frequency among the spots, a q value < 0.05  \cr
 and  a priori probability value < 0.4 was chosen to represent the cluster. \cr
 First run csbl1 after csbl2 in two different R session.
}
\usage{
csbl2(nome_ris, percorso, n_lim = 5, file_set, nome_GO)
}

\arguments{
  \item{nome_ris}{
Filename of output.
}
  \item{percorso}{
Directory where to put the output.
}
  \item{n_lim}{
Minimun number of ids that belongs to a cluster, (default 5).
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
GO_pul_filename \cr
Ids with GO terms, one id per line. \cr \cr
Filename_idclu \cr
Ids divided into clusters, one cluster per line. \cr \cr
Filenames_tab_go \cr
Table with GO terms of selected clusters, frequencies and q value. \cr \cr
Filenames_tab_DEF \cr
Table with GO terms of selected clusters and percentual. \cr \cr
Directory enr \cr
there are files of enrichment analysis done on clusters data. \cr \cr
Bargraph of selected clusters. \cr \cr
Bargraph of selected clusters with table.

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
#id=readLines("/home/marco/biosegen_array/definitivo/setdiff_GRM")
#nome_ris="setdiff"
#percorso="/home/marco/biosegen_array/definitivo/csbl_pro/sim/"
#file_set="/home/marco/funzioni_R/agriGO/prob_tab_sim_slimGRM_def"
#nome_GO="/home/marco/funzioni_R/agriGO/slimGO_GRM_csbl_pul"
#csbl2(nome_ris,percorso,n_lim=5,file_set,nome_GO)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
