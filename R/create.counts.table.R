create.counts.table <- function(corr.summary = NULL) { 

	# establish frequency count table for genes across cancer types
	pancancer.corr.features <- NULL;
	cancer.subtypes.names <- NULL;
	for (cancer.type in names(corr.summary)) {
		pancancer.corr.features <- union(
			pancancer.corr.features, 
			as.vector(unlist(corr.summary[[cancer.type]]$correlated.genes.subtypes))
			);

		# assemble list of subtype names
		for (subtype.name in names(corr.summary[[cancer.type]]$correlated.genes.subtypes)) {
			cancer.subtypes.names <- c(
				cancer.subtypes.names, 
				paste0(cancer.type, ".", subtype.name)
				);
			}
		}
	counts.table <- rep(0, length(pancancer.corr.features));
	names(counts.table) <- pancancer.corr.features;
	counts.table.matrix <- matrix(
		data = 0,
		nrow = length(pancancer.corr.features),
		ncol = length(cancer.subtypes.names),
		dimnames = list(
			pancancer.corr.features,
			cancer.subtypes.names
			)
		);

	for (cancer.type in names(corr.summary)) {
		for (subtype.name in names(corr.summary[[cancer.type]]$correlated.genes.subtypes)) {
			correlated.features <- corr.summary[[cancer.type]]$correlated.genes.subtypes[[subtype.name]];
			if (length(correlated.features) > 0) {
				counts.table[correlated.features] <- counts.table[correlated.features] + 1;
				counts.table.matrix[correlated.features, paste0(cancer.type, ".", subtype.name)] <- 1;
				}
			}
		}

	# sort by frequency of occurrence
	counts.table <- sort(counts.table, decreasing = TRUE);

	# add summary column
	counts.table.matrix <- cbind(
		counts.table.matrix,
		"Frequency" = apply(
			counts.table.matrix,
			1,
			sum
			)
		);

	# sort table as well
	counts.table.matrix <- counts.table.matrix[names(counts.table), , drop = FALSE];

	return(	counts.table.matrix );
	}
