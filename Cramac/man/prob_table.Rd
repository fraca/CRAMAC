\name{prob_table}
\alias{prob_table}

\title{
prob_table
}
\description{
Create similarity and enrichment probability table for csbl.go.
}
\usage{
prob_table(nome, tab_in, dir)
}

\arguments{
  \item{nome}{
Filename output.
}
  \item{tab_in}{
Table with ids and GO. \cr
Zm.6307.1.A1_at	GO:0009408 \cr
Zm.6307.1.A1_at	GO:0009651 \cr
Zm.6309.1.A1_at	GO:0019482 \cr
Zm.6309.1.A1_at	GO:0046251 \cr
Zm.631.1.S1_at	GO:0008094 \cr
Zm.6310.1.A1_at	GO:0016787
}
  \item{dir}{
Directory where to place files.
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
File with enrichment probability table. \cr
File with enrichment probability table. \cr
File with ids and GO, one id per line.
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
# tab=read.table("/home/marco/Microarray/csbl/Maize_spp_AgriGO",
#header=TRUE,stringsAsFactors=FALSE)
# tab=tab[which(tab[[1]]=="maizeAffy"),]
# tab_in=tab[,c(2,3)]
# 
# nome="Affy_AgriGO"
# 
# prob_table(nome,tab_in)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
