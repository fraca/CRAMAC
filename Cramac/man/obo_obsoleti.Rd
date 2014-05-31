\name{obo_obsoleti}
\alias{obo_obsoleti}

\title{
obo_obsoleti
}
\description{
Find function for obsolete GO term.
}
\usage{
obo_obsoleti(obo_tot, nome)
}

\arguments{
  \item{obo_tot}{
File obo to scan.
}
  \item{nome}{
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
File with obsolete GO term and valid GO term.
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
#obo_tot="/home/marco/funzioni_R/csblGRM/gene_ontology_ext.obo"
#nome="/home/marco/funzioni_R/csblGRM/obsoleti"
#obo_obsoleti(obo_tot,nome)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
