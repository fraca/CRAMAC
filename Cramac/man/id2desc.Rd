\name{id2desc}
\alias{id2desc}

\title{
id2desc
}
\description{
the function gives description from ids.
}
\usage{
id2desc(id, tipo, fout = "/tmp/desc_out")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
Vector of ids.
}
  \item{tipo}{
Type of ids: \cr
"GRM" maizesequence	\cr
"array" maizearray	\cr
"affy" affymetrix	\cr
}
  \item{fout}{
Filename output file, (default "/tmp/desc_out").
}
}
\details{
Files needed: \cr
array /home/marco/maizeseq/maizearray.v4.annotation.txt \cr
affy /home/marco/maizeseq/Maize.na32.annot.csv
GRM /home/marco/maizeseq/display.txt
GRM /home/marco/maizeseq/xref.txt
}
\value{
Vector of descriptions and file.
}

\author{
Marco Fracassetti
}

\examples{
#id=readLines("/home/marco/Microarray/nuovo_16_02/id_O6locus_B73_5a")
#tipo="GRM"
#ris=id2desc(id,tipo)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
