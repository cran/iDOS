\name{get.top.features}
\alias{get.top.features}
\title{get.top.features}
\description{
Prioritise top features satisfying the criteria specified by various parameters described below
}
\usage{
get.top.features(
  DE.features = NULL, 
  cna.data.fractions = NULL, 
  mRNA.FC.up = 0, 
  mRNA.FC.down = 0, 
  mRNA.p = 0.05, 
  mRNA.top.n = NULL, 
  cna.fractions.gain = 0.2, 
  cna.fractions.loss = 0.2
  )
}
\arguments{
  \item{DE.features}{Matrix containing differentially expressed features with two columns: FC and P. P may contain adjusted P or raw}
  \item{cna.data.fractions}{Feature by cancer type matrix with CNA fractions}
  \item{mRNA.FC.up}{Log2 fold change threshold for selecting over-expressed features}
  \item{mRNA.FC.down}{Log2 fold change threshold for selecting under-expressed features}
  \item{mRNA.p}{P value threshold for selecting significantly differentially expressed features. Mutually exclusive to \code{mRNA.top.n}}
  \item{mRNA.top.n}{Top n differentially expressed features satisfying each of the fold change criteria. Mutually exclusive to \code{mRNA.p}}
  \item{cna.fractions.gain}{Threshold for selecting copy number gain/amplifications}
  \item{cna.fractions.loss}{Threshold for selecting copy number losses}
}
\value{
Vector of top features
}
\author{
Syed Haider
}
\examples{

# load test data
x <- get.test.data(data.types = c("mRNA.T", "mRNA.N", "CNA"));

# list of features to be assessed for differential expression
feature.ids <- rownames(x$mRNA.T$BLCA);

# get differentially expressed features
DE.results <- find.DE.features(
  exp.data.T = x$mRNA.T, 
  exp.data.N = x$mRNA.N, 
  feature.ids = feature.ids,
  test.name = "t.test"
  );

# get top features
top.features <- get.top.features(
  DE.features = cbind("FC" = DE.results[, 1], "P" = DE.results[, 2]),
  cna.data.fractions = x$CNA.fractions$BLCA, 
  mRNA.FC.up = 0.25, 
  mRNA.FC.down = 0.25, 
  mRNA.p = 0.05, 
  mRNA.top.n = NULL, 
  cna.fractions.gain = 0.2, 
  cna.fractions.loss = 0.2
  );

}
\keyword{ candidate drivers }
