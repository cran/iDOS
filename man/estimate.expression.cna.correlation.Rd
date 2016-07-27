\name{estimate.expression.cna.correlation}
\alias{estimate.expression.cna.correlation}
\title{estimate.expression.cna.correlation}
\description{
Estimate subtype specific correlation between mRNA and CNA profiles
}
\usage{
estimate.expression.cna.correlation(
  exp.data = NULL, 
  cna.data.log2 = NULL, 
  corr.threshold = 0.3, 
  corr.direction = "two.sided", 
  subtypes.metadata = NULL, 
  feature.ids = NULL, 
  cancer.type = NULL, 
  data.dir = NULL, 
  graphs.dir = NULL
  )
}
\arguments{
  \item{exp.data}{Feature by sample mRNA abundance matrix}
  \item{cna.data.log2}{Feature by sample CNA log ratio matrix}
  \item{corr.threshold}{Threshold for Spearman's Rho to consider a feature as candidate driver}
  \item{corr.direction}{Whether to include positively (greater), negatively (less) or both (two.sided) correlated features. Defaults to \code{two.sided}}
  \item{subtypes.metadata}{Subtypes metadata list of lists. Must contain at least one subtype specific samples using list \code{subtype.samples.list}. If no subtypes are present, specify list element "All" with all samples}
  \item{feature.ids}{Vector of features to be used to estimate correlation}
  \item{cancer.type}{Name of the cancer type or dataset}
  \item{data.dir}{Path to output directory where mRNA and CNA correlation statistics will be stored}
  \item{graphs.dir}{Path to graphs directory}
}
\value{
A list of lists containing correlated features per cancer subtype
}
\author{
Syed Haider
}
\examples{

# load test data
x <- get.test.data(data.types = c("mRNA.T", "CNA"));

# temporary output directory
tmp.output.dir <- tempdir();

# estimate mRNA and CNA correlation
correlated.features <- estimate.expression.cna.correlation(
  exp.data = x$mRNA.T$BLCA, 
  cna.data.log2 = x$CNA.log2$BLCA, 
  corr.threshold = 0.3, 
  corr.direction = "two.sided", 
  subtypes.metadata = list(
    "subtype.samples.list" = list("All" = colnames(x$mRNA.T$BLCA))
    ), 
  feature.ids = rownames(x$mRNA.T$BLCA), 
  cancer.type = "BLCA", 
  data.dir = paste(tmp.output.dir, "/data/BLCA/", sep = ""),
  graphs.dir = paste(tmp.output.dir, "/graphs/BLCA/", sep = "")
  );

}
\keyword{ correlation }
