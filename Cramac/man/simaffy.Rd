\name{simaffy}
\alias{simaffy}
\title{
simaffy
}
\description{
The function analyzes affymetrix data.
}
\usage{
simaffy(dir, conf, fold_change = 0.58, p_value = 0.05, 
nome_covdesc = "covdesc.txt", heat_flag = FALSE, graph_flag = FALSE)
}

\arguments{
  \item{dir}{
 Directory where the CEL files are placed.
}
  \item{conf}{
 List of the comparison that will be done.
}
  \item{fold_change}{
 Fold change for significant spot, (default 0.58).
}
  \item{p_value}{
 p value for significant spot, (default 0.05).
}
  \item{nome_covdesc}{
 Name of the file covdesc, (default "covdesc.txt").
}
  \item{heat_flag}{
  If TRUE heatmap will be done, (default FALSE).
}
  \item{graph_flag}{
  If TRUE heatmap will be done, (default FALSE).
}
}
\details{
  Ue the package simpleaffy.
}
\value{
  Table of comparisons. \cr
  Pre-processing image if graph_flag = TRUE. \cr
  heatmap if heat_flag = TRUE
}

\author{
Marco Fracassetti
}

\examples{
# dir="/home/marco/fusarium/array/"
# conf=list(c("fus","open"),c("fus","h2o"),c("h2o","open"))
# simaffy(dir,conf)
# fus_open=read.table("/home/marco/fusarium/array/ris_fus_open.txt")
# fus_open=rownames(fus_open)
# fus_h2o=read.table("/home/marco/fusarium/array/ris_fus_h2o.txt")
# fus_h2o=rownames(fus_h2o)
# h2o_open=read.table("/home/marco/fusarium/array/ris_h2o_open.txt")
# h2o_open=rownames(h2o_open)
# ris=setdiff(intersect(fus_open,fus_h2o),h2o_open)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }