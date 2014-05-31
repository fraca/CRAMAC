\name{GRM2array}
\alias{GRM2array}

\title{
GRM2array
}
\description{
From maizesequence id (GRM) to array id.
}
\usage{
GRM2array(id_maize, tab_conf, score_length = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id_maize}{
 Vector of GRM id.
}
  \item{tab_conf}{
 File with ids table. \cr
ID_affy	ID_maizeseq	score_length	pvalue \cr
Zm.1000.1.A1_at	AC231180.2_FGT006	2.85233918128655	0 \cr
Zm.10015.1.A1_at	GRMZM2G109680_T01	1.81548599670511	0 \cr
Zm.10015.1.A1_at	GRMZM2G109680_T04	1.0164744645799	2e-175 \cr
Zm.10009.1.A1_at	GRMZM2G133629_T01	1.50325097529259	0 \cr
Zm.10009.1.A1_at	GRMZM2G133629_T02	1.43302990897269	0 \cr
Zm.10003.1.A1_at	GRMZM2G402977_T01	1.75844594594595	0
}
  \item{score_length}{
 Minimum score blast in tab_conf, (default NULL).
}
}
\details{
 The tab_conf have be created by blast the array spots against maizesequence transcript database. \cr
 Only results with score greater than 1 and p value less 0.05 have been taken into account. \cr
 created tables: \cr
 /home/marco/maizeseq/affy_B73_5a_cdna_pul \cr
 /home/marco/maizeseq/affy_B73_5b_cdna_pul \cr
 /home/marco/maizeseq/maizearray_B73_5a_cdna_pul \cr
 /home/marco/maizeseq/maizearray_B73_5b_cdna_pul \cr
}
\value{
Table of GRM ids and array ids.
}

\author{
Marco Fracassetti
}

\examples{
# tab_conf="/home/marco/maizeseq/maizearray_B73_5b_cdna_pul"
# id_maize=c("GRMZM2G370852","GRMZM2G404855")
# tab_ris=GRM2array(id_maize,tab_conf)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
