find.DE.features <- function(exp.data.T = NULL, exp.data.N = NULL, feature.ids = NULL, test.name = "t.test") {

	# TODO:
	# sanity checks to go here:

	# set the output columns for various test.name
	test.output.columns <- list(
		"t.test" = c(".FC", ".Padj"),
		"wilcox.test" = c(".FC", ".Padj"),
		"var.test" = c(".FC", ".sigma.T", ".sigma.N", ".Padj")
		);

	all.features <- NULL;
	for (cancer.type in names(exp.data.T)) {
		all.features <- union(all.features, rownames(exp.data.T[[cancer.type]]));
		}

	# limit to features which are present in at least one study
	feature.ids <- intersect(feature.ids, all.features);

	# t-test find cancer-specific differentially expressed genes
	test.stats <- matrix(
		data = NA, 
		nrow = length(feature.ids), 
		ncol = length(exp.data.T) * length(test.output.columns[[test.name]]),
		dimnames = list(
			feature.ids,
			paste(
				rep(
					names(exp.data.T), 
					each = length(test.output.columns[[test.name]])
					), 
				test.output.columns[[test.name]], 
				sep = ""
				)
			)
		);

	for (cancer.type in names(exp.data.T)) {

		# common features between sabatini and this.cancer type
		common.features <- intersect(
			feature.ids, 
			intersect(
				rownames(exp.data.T[[cancer.type]]),
				rownames(exp.data.N[[cancer.type]])
				)
			);

		cat("\nCancer: ", cancer.type, "\tNo. of common features available for differential analysis:\t", length(common.features));

		# make column names
		cancer.type.FC <- paste(cancer.type, ".FC", sep = "");
		cancer.type.P <- paste(cancer.type, ".Padj", sep = "");
		cancer.type.sd.T <- paste(cancer.type, ".sigma.T", sep = "");
		cancer.type.sd.N <- paste(cancer.type, ".sigma.N", sep = "");

		for (feature in common.features) {

			# store fold change
			test.stats[feature, cancer.type.FC] <- round(
				mean(exp.data.T[[cancer.type]][feature, ]) - mean(exp.data.N[[cancer.type]][feature, ]), 
				digits = 3
				);

			# run statistical test and store P value
			if (test.name == "t.test") {
				test.stats[feature, cancer.type.P] <- t.test(
					x = exp.data.T[[cancer.type]][feature, ],
					y = exp.data.N[[cancer.type]][feature, ]
					)$p.value;
				}
			else if(test.name == "wilcox.test") {
				test.stats[feature, cancer.type.P] <- wilcox.test(
					x = exp.data.T[[cancer.type]][feature, ],
					y = exp.data.N[[cancer.type]][feature, ]
					)$p.value;
				}
			else if(test.name == "var.test") {
				test.stats[feature, cancer.type.P] <- var.test(
					x = exp.data.T[[cancer.type]][feature, ],
					y = exp.data.N[[cancer.type]][feature, ]
					)$p.value;
				test.stats[feature, cancer.type.sd.T] <- sd(exp.data.T[[cancer.type]][feature, ], na.rm = T);
				test.stats[feature, cancer.type.sd.N] <- sd(exp.data.N[[cancer.type]][feature, ], na.rm = T);
				}
			else {
				stop("\nDieing gracefully bcoz test.name not recognised");
				}
			}

		# adjust P values
		test.stats[, cancer.type.P] <- p.adjust(
			test.stats[, cancer.type.P], 
			method = "BH"
			);
		}

	return (test.stats);
	}
