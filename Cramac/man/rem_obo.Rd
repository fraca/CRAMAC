\name{rem_obo}
\alias{rem_obo}

\title{
rem_obo
}
\description{
Removing function for obsolete GO term.
}
\usage{
rem_obo(file_in, obs, file_out)
}

\arguments{
  \item{file_in}{
Filename of table with ids and GO. \cr
Zm.6307.1.A1_at	GO:0009408 \cr
Zm.6307.1.A1_at	GO:0009651 \cr
Zm.6309.1.A1_at	GO:0019482 \cr
Zm.6309.1.A1_at	GO:0046251 \cr
Zm.631.1.S1_at	GO:0008094 \cr
Zm.6310.1.A1_at	GO:0016787
}
  \item{obs}{
File whit obsolete GO term, (created with obo_obsoleti).
}
  \item{file_out}{
Filename output.
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
File with non obsolete GO terms, ready for prob_table.
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
#file_in="/home/marco/funzioni_R/csbl_old/csblGRM/GO_GRM"
#obs="/home/marco/funzioni_R/csbl_old/csblGRM/obsoleti"
#file_out="/home/marco/funzioni_R/csbl_old/csblGRM/GO_GRM_pul"
#rem_obo(file_in, obs, file_out)


}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
