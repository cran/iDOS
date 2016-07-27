load.datasets <- function(data.dir = "./", metadata = NULL, data.types = c("mRNA.T", "ann")) {

	# TODO:
	# sanity checks to go here:

	ret.data <- list();

	for (cancer.type in rownames(metadata)) {

		# read mRNA data
		if ("mRNA.T" %in% data.types) {

			# read both Tumour dataset
			ret.data[["mRNA.T"]][[cancer.type]] <- as.matrix(
				read.table(
					file = paste(data.dir, cancer.type, "/tumour/", metadata[cancer.type, "mRNA.T"], sep = ""),
					row.names = 1, 
					header = T, 
					sep = "\t", 
					stringsAsFactors = FALSE
					)
				);

			# remove genes whereby > 75% samples have zero reads
			ret.data[["mRNA.T"]][[cancer.type]] <- ret.data[["mRNA.T"]][[cancer.type]][
				apply(
					ret.data[["mRNA.T"]][[cancer.type]],
					1, 
					FUN = function(x) {
						if (length(which(x == 0)) > (length(x)*0.75)) {
							return (FALSE);
							}
						else {
							return (TRUE);
							}
						}
					), 
				];
			}

		if ("mRNA.N" %in% data.types) {

			# read both Normal dataset

			ret.data[["mRNA.N"]][[cancer.type]] <- as.matrix(
				read.table(
					file = paste(data.dir, cancer.type, "/normal/", metadata[cancer.type, "mRNA.N"], sep = ""),
					row.names = 1, 
					header = T, 
					sep = "\t", 
					stringsAsFactors = FALSE
					)
				);

			# remove genes whereby > 75% samples have zero reads
			ret.data[["mRNA.N"]][[cancer.type]] <- ret.data[["mRNA.N"]][[cancer.type]][
				apply(
					ret.data[["mRNA.N"]][[cancer.type]],
					1, 
					FUN = function(x) {
						if (length(which(x == 0)) > (length(x)*0.75)) {
							return (FALSE);
							}
						else {
							return (TRUE);
							}
						}
					), 
				];
			}


		# read CNA log2 and calls
		if ("CNA" %in% data.types) {

			ret.data[["CNA.log2"]][[cancer.type]] <- as.matrix(
				read.table(
					file = paste(data.dir, cancer.type, "/tumour/", metadata[cancer.type, "CNA.log2"], sep = ""),
					row.names = 1, 
					header = T, 
					sep = "\t", 
					stringsAsFactors = FALSE
					)
				);

			ret.data[["CNA.calls"]][[cancer.type]] <- as.matrix(
				read.table(
					file = paste(data.dir, cancer.type, "/tumour/", metadata[cancer.type, "CNA.calls"], sep = ""),
					row.names = 1, 
					header = T, 
					sep = "\t", 
					stringsAsFactors = FALSE
					)
				);

			ret.data[["CNA.fractions"]][[cancer.type]] <- read.table(
				file = paste(data.dir, cancer.type, "/tumour/", metadata[cancer.type, "CNA.fractions"], sep = ""),
				row.names = 1, 
				header = T, 
				sep = "\t", 
				stringsAsFactors = FALSE
				);
			}

		# read mutations data
		if ("mutations" %in% data.types) {

			ret.data[["mutations"]][[cancer.type]] <- as.matrix(
				read.table(
					file = paste(data.dir, cancer.type, "/tumour/", metadata[cancer.type, "mutations"], sep = ""),
					row.names = 1, 
					header = T, 
					sep = "\t", 
					stringsAsFactors = FALSE
					)
				);
			}

		# read annotation data
		if ("ann" %in% data.types) {

			ret.data[["ann"]][[cancer.type]] <- read.table(
				file = paste(data.dir, cancer.type, "/tumour/", metadata[cancer.type, "annotations"], sep = ""),
				row.names = 1, 
				header = T, 
				sep = "\t", 
				stringsAsFactors = FALSE
				);
			}

		# limit to common samples
		intersectable.data.types <- setdiff(names(ret.data), c("mRNA.N", "CNA.fractions"));
		if (length(intersectable.data.types) > 1) {

			# initialise to samples of first datatype
			common.samples <- colnames(
				ret.data[[
					setdiff(names(ret.data), c("mRNA.N", "CNA.fractions", "ann"))[1]
					]][[cancer.type]]
				);

			# find common samples
			for (data.type in names(ret.data)) {
				if (data.type %in% c("mRNA.N", "CNA.fractions")) { 
					next; 
					}
				else if (data.type == "ann") { 
					common.samples <- intersect(common.samples, rownames(ret.data[[data.type]][[cancer.type]]));
					}
				else {
					common.samples <- intersect(common.samples, colnames(ret.data[[data.type]][[cancer.type]]));
					}
				}

			# limit to common samples
			for (data.type in names(ret.data)) {
				if (data.type %in% c("mRNA.N", "CNA.fractions")) {
					next;
					}
				else if (data.type == "ann") {
					ret.data[[data.type]][[cancer.type]] <- ret.data[[data.type]][[cancer.type]][common.samples, ];
					}
				else {
					ret.data[[data.type]][[cancer.type]] <- ret.data[[data.type]][[cancer.type]][, common.samples];
					}
				}
			}
		}

	return (ret.data);
	}
