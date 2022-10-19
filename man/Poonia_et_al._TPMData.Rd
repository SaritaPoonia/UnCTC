\name{Poonia_et_al._TPMData}
\alias{Poonia_et_al._TPMData}
\docType{data}
\title{Poonia et al.TPM CTC Data}
\description{Single cell RNA-seq TPM matrix of Circulating tumor cells
from poonia et al. data}
\usage{data("Poonia_et_al._TPMData")}
\format{
  The format is:
 num [1:23732, 1:72] 0 0 2.74 0 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:23732] "A1BG" "A1BG-AS1" "A1CF" "A2LD1" ...
  ..$ : chr [1:72] "1851009049_CS5_S22" "1851013039_CS34_S53" "1851009049_CS14_S17" "1851009049_CS9_S25" ...
}

\examples{
data(Poonia_et_al._TPMData)
## maybe str(Poonia_et_al._TPMData) ; plot(Poonia_et_al._TPMData) ...
}
\keyword{datasets}
