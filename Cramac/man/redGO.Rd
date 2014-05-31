\name{redGO}
\alias{redGO}

\title{
redGO
}
\description{
The GO dataset will be reduce by GOslim.
}
\usage{
redGO(file_in, nome, file_obo, onto = "")
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
  \item{nome}{
Filename output.
}
  \item{file_obo}{
file obo with GOslim terms.
}
  \item{onto}{
Ontologies used: \cr
"MF" molecular function \cr
"CC" cellular compartiment \cr
"BP" biological process \cr
"" all \cr
(default "")
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
File of table with ids and GOslim term.
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

# file_in="/home/marco/funzioni_R/csblGRM/GO_GRM_pul_BP"
# nome="/home/marco/funzioni_R/csblGRM/slimGO_GRM_pul_BP"
# file_obo="/home/marco/funzioni_R/csblGRM/map_slimplant"
# onto="BP"
#  redGO(file_in,nome,onto="BP")
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
