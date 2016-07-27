get.top.features <- function(DE.features = NULL, cna.data.fractions = NULL, mRNA.FC.up = 0, mRNA.FC.down = 0, mRNA.p = 0.05, mRNA.top.n = NULL, cna.fractions.gain = 0.2, cna.fractions.loss = 0.2) {

	# sanity checks
	if (!is.null(mRNA.p) && !is.null(mRNA.top.n)) {
		stop("\nDieing gracefully, cannot use both mRNA.p and mRNA.top.n. set one to NULL");
		}

	slimmed.features <- NULL;

	# apply mRNA.p based filter
	if (!is.null(mRNA.p)) {

		# get over-expressed features
		slimmed.features <- rownames(DE.features)[
			which(
				DE.features[, "FC"] > mRNA.FC.up
				&
				DE.features[, "P"] < mRNA.p
				)
			];
	
		# get under-expressed features
		slimmed.features <- c(
			slimmed.features,
			rownames(DE.features)[
				which(
					DE.features[, "FC"] < mRNA.FC.down
					&
					DE.features[, "P"] < mRNA.p
					)
				]
			);
		}
	# apply mRNA.top.n based filter
	else if(!is.null(mRNA.top.n)) {

		if (mRNA.top.n <= 0) {
			stop("\nDieing gracefully, mRNA.top.n must be > 0");
			}

		# get top over-expressed features
		DE.features <- DE.features[order(DE.features[, "FC"], decreasing = TRUE), ];
		slimmed.features <- rownames(DE.features[1:mRNA.top.n, ])[
			which(DE.features[1:mRNA.top.n, "FC"] > mRNA.FC.up)
			];

		# get top under-expressed features
		DE.features <- DE.features[order(DE.features[, "FC"], decreasing = FALSE), ];
		slimmed.features <- c(
			slimmed.features,
			rownames(DE.features[1:mRNA.top.n, ])[
				which(DE.features[1:mRNA.top.n, "FC"] < mRNA.FC.down)
				]
			);
		}
	else {
		}

	# limit to genes with fraction gained
	which.gain <- which(cna.data.fractions[slimmed.features, "All_Gain"] > cna.fractions.gain);

	# limit to genes with fraction loss
	which.loss <- which(cna.data.fractions[slimmed.features, "All_Loss"] > cna.fractions.loss);

	if (length(c(which.gain, which.loss)) > 0) {
		slimmed.features <- unique(slimmed.features[c(which.gain, which.loss)]);
		}
	else {
		slimmed.features <- NULL;
		}

	return (slimmed.features);
	}