estimate.expression.cna.correlation <- function(exp.data = NULL, cna.data.log2 = NULL, corr.threshold = 0.3, corr.direction = "two.sided", subtypes.metadata = NULL, feature.ids = NULL, cancer.type = NULL, data.dir = NULL, graphs.dir = NULL) {

	# sanity test
	if (!identical(colnames(exp.data), colnames(cna.data.log2))) {
		stop("\nDieing gracefully bcoz colnames(exp.data) != colnames(cna.data.log2)");
		}

	if (!identical(rownames(exp.data), rownames(cna.data.log2))) {
		stop("\nDieing gracefully bcoz rownames(exp.data) != rownames(cna.data.log2)");
		}

	# create output directory
	if (!file.exists(data.dir)) {
		dir.create(data.dir, recursive = TRUE);
		}

	if (!file.exists(graphs.dir)) {
		dir.create(graphs.dir, recursive = TRUE);
		}

	# whether to skip venn, in case of empty sets
	plot.venn <- TRUE;

	# get subtype specific samples
	subtype.samples.list <- subtypes.metadata[["subtype.samples.list"]];

	# estimate correlation between mRNA and copy number
	corr.data <- matrix(
		data = NA, nrow = length(feature.ids), ncol = 3, 
		dimnames= list(feature.ids, c("rho", "P", "Q"))
		);
	corr.data.subtypes <- list();
	corr.threshold.genes <- vector();
	corr.threshold.genes.subtypes <- list();
	for (subtype.name in names(subtype.samples.list)) {
		corr.threshold.genes.subtypes[[subtype.name]] <- vector();
		}

	if (length(feature.ids) > 0) {

		# traverse over all subtypes
		for (subtype.name in names(subtype.samples.list)) {

			cat("\n[Correlation] mRNA v CNA: ", cancer.type, subtype.name);

			corr.data.subtype <- corr.data;

			for (gene.name in feature.ids) {

				# cat("\n\t processing gene:", gene.name);

				# estimate spearman correlation
				corr.tmp <- cor.test(
					x = exp.data[gene.name, subtype.samples.list[[subtype.name]]], 
					y = cna.data.log2[gene.name, subtype.samples.list[[subtype.name]]], 
					method = "spearman"
					);

				corr.data.subtype[gene.name, "rho"] <- corr.tmp$estimate;
				corr.data.subtype[gene.name, "P"] <- corr.tmp$p.value;
				}

			# round off to 3 digits
			corr.data.subtype[, "rho"] <- round(corr.data.subtype[, "rho"], digits = 3);
		
			# adjust P values
			corr.data.subtype[, "Q"] <- p.adjust(corr.data.subtype[, "P"], method = "BH");

			# store subtype-specific correlations in list
			corr.data.subtypes[[subtype.name]] <- corr.data.subtype;

			# correlation threshold genes
			if (corr.direction == "two.sided") {
				which.genes <- which(abs(corr.data.subtype[, "rho"]) > corr.threshold);
				}
			else if (corr.direction == "greater") {
				which.genes <- which(corr.data.subtype[, "rho"] > corr.threshold);
				}
			else if (corr.direction == "less") {
				which.genes <- which(corr.data.subtype[, "rho"] < corr.threshold);
				}
			else {
				stop("\nDieing gracefully bcoz corr.direction is invalid");
				}

			if (length(which.genes) > 0) { 
				corr.threshold.genes <- union(
					corr.threshold.genes, 
					rownames(corr.data.subtype)[which.genes]
					);
				corr.threshold.genes.subtypes[[subtype.name]] <- rownames(corr.data.subtype)[which.genes];
				}
			else {
				plot.venn <- FALSE;
				corr.threshold.genes.subtypes[[subtype.name]] <- NULL;
				cat("\n\tvenn diagram will be not be plotted as empty correlation set for this subtype");
				}

			# write to file system
			write.table(
				x = corr.data.subtype[order(corr.data.subtype[, "rho"], decreasing = TRUE), , drop = FALSE],
				file = paste(data.dir, "mRNA_abundance_cna_correlation__", subtype.name, ".txt", sep = ""),
				row.names = TRUE,
				col.names = NA,
				sep = "\t",
				quote = FALSE
				);
			}

		# limit to highly correlated genes only, and if there are at least 2 genes
		if (length(corr.threshold.genes) > 1 && cancer.type == "Metabric" && plot.venn) {

			# make a 5-way venn of PAM50 subtypes
			venn.diagram(
				x = list(
					"Normal-like" = corr.threshold.genes.subtypes[["Normal"]],
					"LuminalA" = corr.threshold.genes.subtypes[["LumA"]],
					"LuminalB" = corr.threshold.genes.subtypes[["LumB"]],
					"Basal" = corr.threshold.genes.subtypes[["Basal"]],
					"Her2" = corr.threshold.genes.subtypes[["Her2"]]
					),
				imagetype = "png",
				filename = paste(graphs.dir, "mRNA_abundance_cna_correlation__venn_PAM50.png", sep = ""),
				col = "black",
				fill = c("forestgreen", "dodgerblue3", "lightskyblue2", "red", "pink"),
				alpha = 0.50,
				fontface = "bold",
				cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
				 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 1, 1, 1, 1, 1.5),
				cat.col = "black",
				cat.cex = 1.5,
				cat.fontface = "bold",
				margin = 0.23,
				cat.dist = 0.32
				);	
			}
		}

	return (
		list(
			"corr.threshold.genes" = corr.threshold.genes,
			"correlated.genes.subtypes" = corr.threshold.genes.subtypes
			)
		);
	}