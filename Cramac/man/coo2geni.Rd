\name{coo2geni}
\alias{coo2geni}

\title{
coo2geni
}
\description{
 The function give the transcripts present in a genomic region.
}
\usage{
coo2geni(coo, db, fout = "/tmp/out_coo")
}

\arguments{
  \item{coo}{
 Coordinates "chr:start-end".
}
  \item{db}{
 Database to use.
}
  \item{fout}{
 Filename of the transcript sequences file, (default "/tmp/out_coo").
}
}
\details{
 Working set database: \cr
/home/marco/maizeseq/B73_5a_cdna/B73_5a_cdna \cr
Filtered set database: \cr
/home/marco/maizeseq/B73_5b_cdna/B73_5b_cdna \cr
}
\value{
 Table with transcripts and coordinates. \cr
 File with transcript sequences \cr
}

\author{
 Marco Fracassetti
}

\examples{
# coo="8:104,811,553-117,198,140" # pro1 +- 5 centimorgan
# coo="8:110811553-111198140" 
# coo="8:110709532-111356822" #contiene anche tip4.2
# db="/home/marco/maizeseq/B73_5b_cdna/B73_5b_cdna"
# fout="/home/marco/Microarray/pro1_l/cdna_pro1_B73_5_b.fasta"
# db="/home/marco/maizeseq/B73_5a_cdna/B73_5a_cdna"
# fout="/tmp/proooo"
# tab=coo2geni(coo,db,fout) #####tabella con coordinate
# write.table(tab,quote=FALSE,sep="\t",file="/home/marco/Microarray/GRM_locus_O6_B73_5a")
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
