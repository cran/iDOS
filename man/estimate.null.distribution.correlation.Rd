\name{estimate.null.distribution.correlation}
\alias{estimate.null.distribution.correlation}
\title{estimate.null.distribution.correlation}
\description{
Function to estimate probability of observing correlations as high as observed using a feature list of interest
}
\usage{
estimate.null.distribution.correlation(
  exp.data = NULL, 
  cna.data.log2 = NULL, 
  corr.threshold = 0.3, 
  corr.direction = "two.sided", 
  subtypes.metadata = NULL, 
  feature.ids = NULL, 
  observed.correlated.features = NULL, 
  iterations = 50, 
  cancer.type = NULL, 
  data.dir = NULL
  )
}
\arguments{
  \item{exp.data}{Feature by sample mRNA abundance matrix}
  \item{cna.data.log2}{Feature by sample CNA log ratio matrix}
  \item{corr.threshold}{Threshold for Spearman's Rho to consider a feature as candidate driver}
  \item{corr.direction}{Whether to include positively (greater), negatively (less) or both (two.sided) correlated features. Defaults to \code{two.sided}}
  \item{subtypes.metadata}{Subtypes metadata list. Contains at least subtype specific samples}
  \item{feature.ids}{Vector of features to be used to estimate correlation}
  \item{observed.correlated.features}{List of features that were found to be correlated for subtypes of a given cancer type}
  \item{iterations}{Number of random permutations for estimating p value}
  \item{cancer.type}{Name of the cancer type or dataset}
  \item{data.dir}{Path to output directory where the randomisation results will be stored}
}
\value{
1 if successful
}
\author{
Syed Haider
}
\seealso{
\code{\link{estimate.expression.cna.correlation}}
}
\examples{

# load test data
x <- get.test.data(data.types = c("mRNA.T", "CNA"));

# temporary output directory
tmp.output.dir <- tempdir();

# estimate mRNA and CNA correlation for each cancer/disease type
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

# estimate NULL distribution
estimate.null.distribution.correlation(
  exp.data = x$mRNA.T$BLCA,
  cna.data.log2 = x$CNA.log2$BLCA, 
  corr.threshold = 0.3, 
  corr.direction = "two.sided", 
  subtypes.metadata = list(
    "subtype.samples.list" = list("All" = colnames(x$mRNA.T$BLCA))
    ), 
  feature.ids = rownames(x$mRNA.T$BLCA), 
  observed.correlated.features = correlated.features$correlated.genes.subtypes, 
  iterations = 50, 
  cancer.type = "BLCA", 
  data.dir = paste(tmp.output.dir, "/data/BLCA/", sep = "")
  );

}
\keyword{ randomisation }
\keyword{ Null distribution }