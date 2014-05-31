\name{hypep}
\alias{hypep}

\title{
hypep
}
\description{
The function detects and displays various useful metrics about a protein sequence. \cr
Require EMBOSS installed.
}
\usage{
hypep(ids, db)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ids}{
Vector of ids to use.
}
  \item{db}{
Protein database to use.
}
}

\value{
List: \cr
id ids analyzed. \cr
MW molecular weight. \cr
n_AA number of amino acid. \cr
inc_bodies Improbability of expression in inclusion bodies. \cr
tab_AA Table with percentual of amino acid. \cr
Tiny (A+C+G+S+T) \cr
Small(A+B+C+D+G+N+P+S+T+V) \cr
Aliphatic (A+I+L+V) \cr
Aromatic (F+H+W+Y) \cr
Non-polar (A+C+F+G+I+L+M+P+V+W+Y) \cr
Polar (D+E+H+K+N+Q+R+S+T+Z) \cr
Charged (B+D+E+H+K+R+Z) \cr
Basic (H+K+R) \cr
Acidic (B+D+E+Z) \cr  \cr
HY_OHM mean of OHM  Hydropathy parameters (Sweet & Eisenberg) of all amino acids. \cr
HY_Doo mean of Kyte & Doolittle hydropathy parameters of all amino acids. \cr

}

\author{
Marco Fracassetti
}

\examples{
# ids=c("GRMZM2G129133_T01", "GRMZM2G162450_T01", "GRMZM2G401342_T01", "GRMZM2G163015_T01",
# "AC204530.4_FGT003", "AC190636.3_FGT005","GRMZM2G439784_T01", "GRMZM2G006234_T01",
# "GRMZM2G463968_T01", "GRMZM2G319307_T01")
# db="/home/marco/maizeseq/B73_5b_prot/B73_5b_prot"
# ris=hypep(ids,db)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
