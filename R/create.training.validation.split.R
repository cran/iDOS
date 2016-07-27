create.training.validation.split <- function(exp.data = NULL, ann.data = NULL, seed.number = 51214) {


	# sanity test
	if (!identical(colnames(exp.data), rownames(ann.data))) {
		stop("\nDieing gracefully bcoz colnames(exp.data) != rownames(ann.data)");
		}

	# length of training set
	# 66:34 split, otherwise 50:50
	if (ncol(exp.data) < 200) {
		training.length <- round(ncol(exp.data)*0.66);
		}
	else {
		training.length <- round(ncol(exp.data)*0.5);	
		}

	# set seed
	set.seed(seed.number);
	samples.T.indices <- sample(1:ncol(exp.data), training.length, replace = FALSE);

	samples.T <- colnames(exp.data)[samples.T.indices];
	samples.V <- setdiff(colnames(exp.data), samples.T);

	return(
		list(
			"exp.T" = exp.data[, samples.T],
			"ann.T" = ann.data[samples.T, ],
			"exp.V" = exp.data[, samples.V],
			"ann.V" = ann.data[samples.V, ]
			)
		);
	}
