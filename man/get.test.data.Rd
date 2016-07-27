\name{get.test.data}
\alias{get.test.data}
\title{get.test.data}
\description{
Function to load test data
}
\usage{
get.test.data(data.types = c("mRNA.T", "ann"))
}
\arguments{
  \item{data.types}{Datatypes to be read Valid datatypes are: mRNA.T, mRNA.N, CNA (includes: log2, calls and fractions), annotations}
}
\value{
List of lists containing datasets and respective molecular profiles as matrices
}
\author{
Syed Haider
}
\examples{

x <- get.test.data(data.types = c("mRNA.T", "mRNA.N", "ann"));

}
\keyword{ datasets }

