\name{create.counts.table}
\alias{create.counts.table}
\title{create.counts.table}
\description{
Summary function to collapse the counts of selected (e.g. correlated) features per cancer type into counts table
}
\usage{
create.counts.table(corr.summary = NULL)
}
\arguments{
  \item{corr.summary}{A list object containing subtype specific selected (e.g. correlated) features. This is the list object returned by \code{estimate.expression.cna.correlation}}
}
\value{
A matrix of cancer type specific counts
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

# go through each cancer type iteratively and perform mRNA-CNA correlation analysis
correlated.features <- list();
for (cancer.type in names(x$mRNA.T)) {

  # estimate mRNA and CNA correlation for each cancer/disease type
  correlated.features[[cancer.type]] <- estimate.expression.cna.correlation(
    exp.data = x$mRNA.T[[cancer.type]], 
    cna.data.log2 = x$CNA.log2[[cancer.type]], 
    corr.threshold = 0.3, 
    corr.direction = "two.sided", 
    subtypes.metadata = list(
      "subtype.samples.list" = list("All" = colnames(x$mRNA.T[[cancer.type]]))
	  ), 
    feature.ids = rownames(x$mRNA.T[[cancer.type]]), 
    cancer.type = cancer.type, 
    data.dir = paste(tmp.output.dir, "/data/", cancer.type, sep = ""),
    graphs.dir = paste(tmp.output.dir, "/graphs/", cancer.type, sep = "")
    );
  }

# create counts table across cancer types
counts.table <- create.counts.table(corr.summary = correlated.features);

}
\keyword{ output }
