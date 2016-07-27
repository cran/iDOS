\name{load.datasets}
\alias{load.datasets}
\title{load.datasets}
\description{
Function to load and systemise molecular datasets
}
\usage{
load.datasets(
  data.dir = "./", 
  metadata = NULL, 
  data.types = c("mRNA.T", "ann")
  )
}
\arguments{
  \item{data.dir}{Path to base data directory or directory containing molecular profiles}
  \item{metadata}{Dataset by profile metadata matrix containing file names of the molecular profiles for different datasets}
  \item{data.types}{Datatypes to be read Valid datatypes are: mRNA.T, mRNA.N, CNA (includes: log2, calls and fractions), annotations}
}
\value{
List of lists containing datasets and respective molecular profiles as matrices
}
\author{
Syed Haider
}
\examples{

# locate test data directory which comes with the package
data.dir <- paste(system.file("programdata/testdata/", package = "iDOS"), "/", sep = "");

# read meta data file
metadata <- read.table(
  file = paste(data.dir, "metadata.txt", sep = ""), 
  row.names = 1, 
  header = TRUE, 
  sep = "\t",
  stringsAsFactors = FALSE
  );

x <- load.datasets(
  data.dir = data.dir,
  metadata = metadata,
  data.types = c("mRNA.T", "mRNA.N", "ann")
  );

}
\keyword{ datasets }

