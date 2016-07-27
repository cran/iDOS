get.test.data <- function(data.types = c("mRNA.T", "ann")) {

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
		data.types = data.types
		);

	return (x);

	}