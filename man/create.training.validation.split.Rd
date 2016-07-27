\name{create.training.validation.split}
\alias{create.training.validation.split}
\title{create.training.validation.split}
\description{
Utility function to create random partitions of a dataset into training and validation sets. If samples are < 200, 66:34; otherwise 50:50 partitions are generated between training and validation sets respectively
}
\usage{
create.training.validation.split(
  exp.data = NULL, ann.data = NULL, seed.number = 51214
  )
}
\arguments{
  \item{exp.data}{Feature by sample mRNA abundance matrix}
  \item{ann.data}{Sample by clinical attribute matrix}
  \item{seed.number}{Random seed for sampling}
}
\value{
A list of four matrices expression and two associated clinical matrices (exp.T, ann.T, exp.V and ann.V). One set for training and one for validation
}
\author{
Syed Haider
}
\examples{

# load test data
x <- get.test.data(data.types = c("mRNA.T", "ann"));

# create training and validation sets
partitioned.datasets <- create.training.validation.split(
  exp.data = x$mRNA.T$BLCA, 
  ann.data = x$ann$BLCA, 
  seed.number = 51214
  );

}
\keyword{ output }
