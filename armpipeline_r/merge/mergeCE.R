mergeCE <- function(C, E, m = 2, n = 1, dropPar = FALSE) {
	# This function is the first block of the ARM pipeline. 
	# It expands the matrix of cases (C) adding exposures from the matrix (E) according to a common parameter (for example, a column "date").

	# C in input should have 2 columns, the first one being cases/controls.
	# E can have an arbitrary number of columns. The column of the common parameter that the user wishes to use for merging must contain unique values.

	# IMPORTANT: the class of vectors in the database are important for the correct functioning of the function mergeCE. To be sure to preserve the class when importing the data from a csv file, use the command read.csv("my_exposures.csv",stringsAsFactors=FALSE). 
	
	# If m and/or n have been passed as column names instead of index, the index is found by the function
	if (is.character(m)) {
		m <- match(m, colnames(C))
		if (is.na(m)) 
			stop("parameter tag for cases not found.")
	}

	if (is.character(n)) {
		n <- match(n, colnames(E))
		if (is.na(n)) 
			stop("parameter tag for exposures not found.")
	}

	# Test for uniquness of values in E[,n]
	u <- length(unique(E[, n]))
	if (u < nrow(E)) {
		stop("parameter column in exposure matrix contains non-unique values.")
	}

	# Test that the user selected parameter columns are of the same type
	if (class(E[, n]) != class(C[, m])) {
		stop("selected parameter columns are of different data type.")
	}


	# Store original size of C
	colC <- ncol(C)

	for (i in 1:ncol(E)) {
		labels <- colnames(E)
		type <- class(E[, i])
		C[, labels[i]] <- vector(type, nrow(C))
	}

	flag <- TRUE
	sink("missing_match_report.txt")
	cat("Missing match report\n\n")
	cat("Tag \t Row\n")

	for (i in 1:nrow(C)) {
		ind = match(C[i, m], E[, n])
		if (is.na(ind)) {
			if (flag) 
				warning("Match not found.")
			cat(sprintf('%s \t %s\n', as.character(C[i, m]), as.character(i)))
		} else C[i, (colC + 1):ncol(C)] <- E[ind, ]
	}

	dropind <- colC + n
	C[, dropind] <- NULL
	if (dropPar) 
		C[, m] <- NULL

	sink()
	closeAllConnections()
	return(C)
}


# SCRAP VERSION SAVED FOR THE RECORD
mergeCEscrap <- function(C, E, m = 2, n = 3) {
	# This function is the first block of the ARM pipeline. 
	# It merges the matrix of cases (C) with the matrix of the exposures (E) according to a common parameter (for example, a column "date") and outputs the resulting matrix (O).

	# C in input should have 2 columns, the first one being cases/controls.
	# E can have an arbitrary number of columns. The colomn of the common parameter that the user wishes to use for merging must contain unique values.

	# Create output dataset O (empty) of size C_rows * C+E_cols
	O <- data.frame(matrix(NA, nrow(C), ncol(C) + ncol(E)))

	# Match parameter in C[,m] with E[,n] and copy in O
	for (i in 1:nrow(C)) {
		ind = match(C[i, m], E[, n])
		O[i, ] <- c(C[i, ], E[ind, ])
	}
	return(O)
}

