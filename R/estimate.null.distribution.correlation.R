estimate.null.distribution.correlation <- function(exp.data = NULL, cna.data.log2 = NULL, corr.threshold = 0.3, corr.direction = "two.sided", subtypes.metadata = NULL, feature.ids = NULL, observed.correlated.features = NULL, iterations = 50, cancer.type = NULL, data.dir = NULL) {

	# limit to common features
	exp.data <- exp.data[intersect(rownames(exp.data), rownames(cna.data.log2)), ];
	cna.data.log2 <- cna.data.log2[intersect(rownames(exp.data), rownames(cna.data.log2)), ];

	# sanity test
	if (!identical(colnames(exp.data), colnames(cna.data.log2))) {
		stop("\nDieing gracefully bcoz colnames(exp.data) != colnames(cna.data.log2)");
		}

	if (!identical(rownames(exp.data), rownames(cna.data.log2))) {
		stop("\nDieing gracefully bcoz rownames(exp.data) != rownames(cna.data.log2)");
		}

	if (length(feature.ids) == 0) {
		cat("\nWarning: skipping NULL distribution as length(feature.ids) = 0");
		return ();
		}

	sig.size <- length(feature.ids);
	random.corr.p <- list();

	# get subtype specific samples
	subtype.samples.list <- subtypes.metadata[["subtype.samples.list"]];

	# traverse over all subtypes
	for (subtype.name in names(subtype.samples.list)) {

		cat("\n[ESTIMATING NULL DISTRIBUTION - CORRELATION]:", cancer.type, subtype.name);

		corr.observed <- length(observed.correlated.features[[subtype.name]]);
		cat("\tobserved correlation count: ", corr.observed, "/", sig.size);

		# pool of genes
		gene.pool <- rownames(exp.data);
		gt.count <- 0;

		for (i in 1:iterations) {

			random.genes <- sample(gene.pool, sig.size);

			# estimate rho for cna and mRNA
			genes.rho <- unlist(
				lapply(
					random.genes, 
					FUN = function(gene) { 
						cor.test(
							x = exp.data[gene, subtype.samples.list[[subtype.name]]], 
							y = cna.data.log2[gene, subtype.samples.list[[subtype.name]]], 
							method = "spearman"
							)$estimate;
						}
					)
				);

			# count if random event is greater than observed
			if (corr.direction == "two.sided") {
				if (length(which(abs(genes.rho) > corr.threshold)) > corr.observed) {
					gt.count <- gt.count + 1;
					}
				}
			else if (corr.direction == "greater") {
				if (length(which(genes.rho > corr.threshold)) > corr.observed) {
					gt.count <- gt.count + 1;
					}
				}
			else if (corr.direction == "less") {
				if (length(which(genes.rho < corr.threshold)) > corr.observed) {
					gt.count <- gt.count + 1;
					}
				}
			else {
				stop("\nDieing gracefully bcoz corr.direction is invalid");
				}
			}

		# estimate P value
		random.corr.p[[subtype.name]] <- c(
			"Observed" = corr.observed,
			"Total" = sig.size,
			"P" = (gt.count+1)/(iterations+1)
			);
		}

	# make matrix
	random.corr.p <- do.call(rbind, random.corr.p);

	# write to file system
	write.table(
		random.corr.p,
		file = paste(data.dir, "mRNA_abundance_cna_correlation__null_distribution_", iterations, ".txt", sep = ""),
		row.names = TRUE,
		col.names = NA,
		sep = "\t"
		);

	return ();
	}
