\name{t_mel}
\alias{t_mel}

\title{
t_mel
}
\description{
Calculate the melting temperature TM. \cr
formula: \cr
69.3 + (41 X nGC/nTOT) - (650/nTOT)
}
\usage{
t_mel(primer)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{primer}{
String with primer.
}
}

\value{
Print and return TM.
}

\author{
Marco Fracassetti
}

\examples{
#t_mel("TCCTCCGCTTATTGATATGC")
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
