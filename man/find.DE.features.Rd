\name{find.DE.features}
\alias{find.DE.features}
\title{find.DE.features}
\description{
Funtion to identify differentially expressed/variable features between Tumour (T) and Normal (N) profiles
}
\usage{
find.DE.features(
  exp.data.T = NULL, 
  exp.data.N = NULL, 
  feature.ids = NULL,
  test.name = "t.test"
  )
}
\arguments{
  \item{exp.data.T}{Feature by sample mRNA abundance matrix; tumour samples}
  \item{exp.data.N}{Feature by sample mRNA abundance matrix; normal/baseline samples}
  \item{feature.ids}{Vector of features to be used to estimate correlation}
  \item{test.name}{Specify the statistical test name (exactly as it appears in R). Supported tests are \code{t.test}, \code{wilcox.test}, \code{var.test}}
}
\value{
Feature by cancer type matrix of log2 fold change (T vs N) and adjusted P values. P values are estimated through \code{test.name}
}
\author{
Syed Haider
}
\seealso{
\code{\link{t.test}}, 
\code{\link{wilcox.test}}, 
\code{\link{var.test}}
}
\examples{

# load test data
x <- get.test.data(data.types = c("mRNA.T", "mRNA.N"));

# list of features to be assessed for differential expression
feature.ids <- rownames(x$mRNA.T$BLCA);

DE.results <- find.DE.features(
  exp.data.T = x$mRNA.T, 
  exp.data.N = x$mRNA.N, 
  feature.ids = feature.ids,
  test.name = "t.test"
  );

}
\keyword{ differential features }