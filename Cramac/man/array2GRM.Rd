\name{array2GRM}
\alias{array2GRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
array2GRM
}
\description{
From array id  to maizesequence id (GRM).
}
\usage{
array2GRM(id_array, tab_conf, tab_array = NULL, transcript = FALSE,
 score_length = NULL, file_gramene = "/tmp/X_gramene")
}

\arguments{
  \item{id_array}{
 Vector of array id.
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
  \item{tab_array}{
 Table with array data, rownames same id_array, (default NULL).
}
  \item{transcript}{
 If TRUE the function gives transcript, (default FALSE).
}
  \item{score_length}{
 Minimum score blast in tab_conf, (default NULL).
}
  \item{file_gramene}{
 Filename input for gramene, (default "/tmp/X_gramene").
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
If tab_array equal NULL  table with array id and GRM id, \cr
else list table with array id and GRM id, table with GRM id and array data.

}

\author{
Marco Fracassetti
}

\examples{
#esempio
# page=read.table("/home/marco/Microarray/opaque/prova/array/
# risultati_fin/analisi_fin_11_def/X_PAGE",stringsAsFactors=FALSE,header=TRUE)
# rownames(page)=page[,1]
# page=page[unique(union(which(abs(page[,"d1"])>0.58),
# which(abs(page[,"d2"])>0.58))),c("d1","d2")]
# page2=page[which(abs(page[,"d1"])>0.58),"d1"]
# id_array=rownames(page)
# tab_array=page2
# tab_conf="/home/marco/maizeseq/maizearray_B73_5b_cdna_pul"
# transcript=FALSE
# score_length=NULL
# tab=array2GRM(id_array,tab_conf,tab_array=page)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
