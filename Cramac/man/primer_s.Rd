\name{primer_s}
\alias{primer_s}

\title{
primer_s
}
\description{
Find primers in a sequence.
}
\usage{
primer_s(seq, temp = 60, win_min = 19, win_max = 21, n_primer = 5, fin = TRUE)
}

\arguments{
  \item{seq}{
String with the sequence.
}
  \item{temp}{
Melting temperature (default 60).
}
  \item{win_min}{
Minimum length primer, (default 19).
}
  \item{win_max}{
Maximum length primer, (default 21).
}
  \item{n_primer}{
Number of primers displayed.
}
  \item{fin}{
If TRUE the primer have to start and finish with C or G, (default TRUE).
}
}

\value{
Print the primers with TM.
}

\author{
Marco Fracassetti
}

\examples{
#seq_L="ATTTTGTAATCAATTAAGCGAGGCGATGAGCTCGTCCTCTCTCTTGCACGAGCC"
#primer_s(seq_L)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
