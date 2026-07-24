## allelematch R Package
## v3.0.0
## allelematch: Pairwise matching and identification of unique multilocus genotypes
##
## by Paul Galpern
## License: GPL-2
##
## Last update: 23 July 2026
##
## N.B. Renamed from MicroSatMatch as of v1.1
##
## Functions:
## amDataset()             Produces an input dataset object for allelematch routines
## print.amDataset()       Print method for amDataset objects
## amPreCheck()            Pre-screens an amDataset object for severe missing data loads
## amMatrix()              Produce a dissimilarity matrix
## amPairwise()            Pairwise matching of genotypes
## summary.amPairwise()    Summary method for amPairwise objects
## amCluster()             Clustering of genotypes
## summary.amCluster()     Summary method for amCluster objects
## amAlleleFreq()          Produces allele frequencies from an amDataset object
## amUnique()              Identification of unique genotypes
## amUniqueProfile()       Utility to find optimal parameters for amUnique()
## summary.amUnique()      Summary method for amUnique objects
##
## Requires: dynamicTreeCut
##
## Please see R documentation for full description of function parameters


#### amDataset() ####
amDataset <-
  function(multilocusDataset,
           missingCode = "-99",
           indexColumn = NULL,
           metaDataColumn = NULL,
           ignoreColumn = NULL) {
    newDataset <- list()
    class(newDataset) <- "amDataset"

    if (is.null(dim(multilocusDataset)))
      stop("allelematch: multilocusDataset must be a matrix or a data.frame",
           call. = FALSE)

    if ((
      is.character(indexColumn) ||
      is.character(metaDataColumn) || is.character(ignoreColumn)
    ) &&
    is.null(dimnames(multilocusDataset)[[2]])) {
      stop(
        "allelematch: multilocusDataset does not have dimnames()[[2]] set; use integer column indices or set dimnames()",
        call. = FALSE
      )
    }

    if (!is.null(indexColumn)) {
      if (length(indexColumn) > 1)
        stop("allelematch: only one indexColumn permitted.", call. = FALSE)
      if (is.character(indexColumn)) {
        indexColumnWhich <- which(indexColumn == dimnames(multilocusDataset)[[2]])
        if (length(indexColumnWhich) == 0)
          stop("allelematch: indexColumn does not exist in multilocusDataset",
               call. = FALSE)
      } else {
        indexColumnWhich <- indexColumn
        if (length(indexColumn) &&
            indexColumnWhich > ncol(multilocusDataset)) {
          stop("allelematch: indexColumn does not exist in multilocusDataset",
               call. = FALSE)
        }
      }
    } else {
      indexColumnWhich <- 0
    }

    if (!is.null(metaDataColumn)) {
      if (length(metaDataColumn) > 1)
        stop("allelematch: only one metaDataColumn permitted", call. = FALSE)
      if (is.character(metaDataColumn)) {
        metaDataColumnWhich <- which(metaDataColumn == dimnames(multilocusDataset)[[2]])
        if (length(metaDataColumnWhich) == 0)
          stop("allelematch: metaDataColumn does not exist in multilocusDataset",
               call. = FALSE)
      } else {
        metaDataColumnWhich <- metaDataColumn
        if (metaDataColumnWhich > ncol(multilocusDataset))
          stop("allelematch: metaDataColumn does not exist in multilocusDataset",
               call. = FALSE)
      }
    } else {
      metaDataColumnWhich <- 0
    }

    if (!is.null(ignoreColumn)) {
      if (is.character(ignoreColumn)) {
        ignoreColumnWhich <- which(dimnames(multilocusDataset)[[2]] %in% ignoreColumn)
        if (length(ignoreColumnWhich) != length(ignoreColumn))
          stop(
            "allelematch: one or more ignoreColumn does not exist in multilocusDataset",
            call. = FALSE
          )
      } else {
        ignoreColumnWhich <- ignoreColumn
        if (any(ignoreColumnWhich > ncol(multilocusDataset)))
          stop(
            "allelematch: one or more ignoreColumn does not exist in multilocusDataset",
            call. = FALSE
          )
      }
    } else {
      ignoreColumnWhich <- 0
    }

    columnDataset <- colnames(multilocusDataset)
    # Force everything to character BEFORE matrix conversion
    # If it's a data.frame, convert factors to character strings first
    if (is.data.frame(multilocusDataset)) {
      multilocusDataset <- data.frame(lapply(multilocusDataset, as.character),
                                      stringsAsFactors = FALSE)
    }

    # Now conversion to matrix will preserve the actual strings
    multilocusDataset <- as.matrix(multilocusDataset)

    if (anyNA(multilocusDataset)) {
      multilocusDataset[is.na(multilocusDataset)] <- missingCode
      cat("allelematch: NA data converted to ", missingCode, "\n", sep = "")
    }

    newDataset$index <- multilocusDataset[, indexColumnWhich]
    newDataset$metaData <- multilocusDataset[, metaDataColumnWhich]
    keepTheseColumns <- 1:ncol(multilocusDataset) %in% c(indexColumnWhich, metaDataColumnWhich, ignoreColumnWhich)

    if (sum(!keepTheseColumns) < 3)
      stop("allelematch: at least three data columns are required for allelematch",
           call. = FALSE)

    newDataset$multilocus <- multilocusDataset[, !keepTheseColumns]
    newDataset$multilocus[] <- gsub(" ", "", newDataset$multilocus, fixed = TRUE)
    newDataset$index <- gsub(" ", "", newDataset$index)
    columnDataset <- columnDataset[!keepTheseColumns]

    if (length(newDataset$index) == 0) {
      if (nrow(newDataset$multilocus) > 17576)
        stop(
          "allelematch: too many samples for automatic assignment of index; please provide an index column",
          call. = FALSE
        )
      grid <- expand.grid(LETTERS, LETTERS, LETTERS)
      labelRepository <- paste0(grid$Var3, grid$Var2, grid$Var1)
      newDataset$index <- labelRepository[1:nrow(newDataset$multilocus)]
    }
    if (length(newDataset$metaData) == 0)
      newDataset$metaData <- NULL

    if (is.null(columnDataset)) {
      dimnames(newDataset$multilocus)[[2]] <- paste0("loc", 1:ncol(newDataset$multilocus))
    } else {
      dimnames(newDataset$multilocus)[[2]] <- columnDataset
    }

    if (anyDuplicated(newDataset$index) > 0)
      stop("allelematch: index column should contain a unique identifier for each sample",
           call. = FALSE)

    if (is.na(missingCode)) {
      newDataset$multilocus[is.na(newDataset$multilocus)] <- "NA"
      newDataset$missingCode <- "NA"
    } else {
      newDataset$missingCode <- missingCode
    }
    return(newDataset)
  }


#### print.amDataset() ####
print.amDataset <- function(x, ...) {
  cat("allelematch\namDataset object\n\n")

  # 1. Print the index component
  cat("$index\n")
  print(x$index)
  cat("\n")

  # 2. Print the multilocus matrix
  cat("$multilocus\n")

  # To match the old format (matrix with row/col labels),
  # print the matrix directly rather than converting to data.frame
  print(x$multilocus)

  # 3. Print the missingCode component
  cat("\n$missingCode\n")
  print(x$missingCode)
}


#### amMatrix() ####
amMatrix <- function(amDatasetFocal, missingMethod = 2) {
  if (!inherits(amDatasetFocal, "amDataset"))
    stop(
      "allelematch: amDatasetFocal must be an object of class \"amDataset\"; use amDataset() first",
      call. = FALSE
    )
  if (!(missingMethod %in% c(1, 2)))
    stop("allelematch: missingMethod must equal 1 or 2", call. = FALSE)

  numGenotypes <- nrow(amDatasetFocal$multilocus)
  genotypes <- amDatasetFocal$multilocus
  genotypes[genotypes == amDatasetFocal$missingCode] <- NA
  matchMatrix <- matrix(0, nrow = numGenotypes, ncol = numGenotypes)

  for (col in 1:ncol(genotypes)) {
    col_val <- genotypes[, col]
    match_mat <- outer(col_val, col_val, "==")
    match_mat[is.na(match_mat)] <- FALSE
    matchMatrix <- matchMatrix + match_mat
  }

  missingMultiplier <- if (missingMethod == 1)
    0.25
  else
    0.5
  comparisonMissing <- rowSums(is.na(genotypes)) * missingMultiplier

  simMatrix <- matchMatrix +
    matrix(comparisonMissing,
           nrow = numGenotypes,
           ncol = numGenotypes,
           byrow = TRUE) +
    matrix(comparisonMissing,
           nrow = numGenotypes,
           ncol = numGenotypes,
           byrow = FALSE)

  storage.mode(simMatrix) <- "double"
  dissimMatrix <- 1 - (simMatrix / ncol(genotypes))
  dimnames(dissimMatrix) <- list(amDatasetFocal$index, amDatasetFocal$index)
  class(dissimMatrix) <- "amMatrix"
  return(dissimMatrix)
}


#### amPairwise() ####
amPairwise <-
  function(amDatasetFocal,
           amDatasetComparison = amDatasetFocal,
           alleleMismatch = NULL,
           matchThreshold = NULL,
           missingMethod = 2) {
    if ((!inherits(amDatasetFocal, "amDataset")) ||
        (!inherits(amDatasetComparison, "amDataset"))) {
      stop(
        "allelematch: amDatasetFocal and amDatasetComparison must be an object of class \"amDataset\"; use amDataset() first",
        call. = FALSE
      )
    }
    if (!(missingMethod %in% c(1, 2)))
      stop("allelematch: missingMethod must equal 1 or 2", call. = FALSE)
    if (sum(!(c(
      is.null(alleleMismatch), is.null(matchThreshold)
    ))) != 1)
      stop("allelematch: please specify alleleMismatch OR matchThreshold",
           call. = FALSE)
    if (length(c(alleleMismatch, matchThreshold)) > 1)
      stop(
        "allelematch: please provide a single parameter value for alleleMismatch OR matchThreshold",
        call. = FALSE
      )

    if (!is.null(matchThreshold)) {
      if ((matchThreshold < 0) ||
          (matchThreshold > 1))
        stop("allelematch: matchThreshold must be between 0 and 1", call. = FALSE)
      alleleMismatch <- round((1 - matchThreshold) * ncol(amDatasetFocal$multilocus), 2)
    } else if (!is.null(alleleMismatch)) {
      matchThreshold <- 1 - (alleleMismatch / ncol(amDatasetFocal$multilocus))
    }

    focalGenotypes <- amDatasetFocal$multilocus
    comparisonGenotypes <- amDatasetComparison$multilocus
    indexFocal <- amDatasetFocal$index
    indexComparison <- amDatasetComparison$index
    metaDataFocal <- amDatasetFocal$metaData
    metaDataComparison <- amDatasetComparison$metaData
    columnNames <- dimnames(amDatasetFocal$multilocus)[[2]]

    if (ncol(focalGenotypes) != ncol(comparisonGenotypes))
      stop(
        "allelematch: amDatasetFocal and amDatasetComparison must have the same number of columns / loci",
        call. = FALSE
      )
    if (any(
      dimnames(amDatasetFocal$multilocus)[[2]] != dimnames(amDatasetComparison$multilocus)[[2]]
    ))
      stop(
        "allelematch: amDatasetFocal and amDatasetComparison must have columns with the same names",
        call. = FALSE
      )
    if (amDatasetFocal$missingCode != amDatasetComparison$missingCode)
      stop(
        "allelematch: amDatasetFocal and amDatasetComparison must have same missingCode",
        call. = FALSE
      )

    focalGenotypes[focalGenotypes == amDatasetFocal$missingCode] <- NA
    comparisonGenotypes[comparisonGenotypes == amDatasetComparison$missingCode] <- NA

    numFocalGenotypes <- nrow(focalGenotypes)
    numComparisonGenotypes <- nrow(comparisonGenotypes)
    simMatrix <- matrix(, nrow = numFocalGenotypes, ncol = numComparisonGenotypes)
    pairwiseMatches <- vector("list", numFocalGenotypes)

    missingMultiplier <- if (missingMethod == 1)
      0.25
    else
      0.5
    compMissingSums <- rowSums(is.na(comparisonGenotypes)) * missingMultiplier
    focalMissingSums <- rowSums(is.na(focalGenotypes)) * missingMultiplier

    for (i in 1:numFocalGenotypes) {
      focal_row <- focalGenotypes[i, ]
      match_count <- rowSums(comparisonGenotypes == rep(focal_row, each = numComparisonGenotypes),
                             na.rm = TRUE)

      simMatrix[i, ] <- as.double(match_count + compMissingSums + focalMissingSums[i])
      simMatrix[i, ] <- simMatrix[i, ] / ncol(focalGenotypes)
      pairwiseMatchesWhich <- which(simMatrix[i, ] >= matchThreshold)
      pairwiseMatchesScores <- signif(simMatrix[i, pairwiseMatchesWhich], 2)

      pairwiseMatches[[i]]$focal <- list(
        index = indexFocal[i],
        metaData = metaDataFocal[i],
        multilocus = t(focalGenotypes[i, ])
      )
      dimnames(pairwiseMatches[[i]]$focal$multilocus)[[2]] <- columnNames
      pairwiseMatches[[i]]$focal$flags <- matrix(as.integer(is.na(pairwiseMatches[[i]]$focal$multilocus)) + 1, 1, ncol = length(columnNames))
      pairwiseMatches[[i]]$focal$multilocus[is.na(pairwiseMatches[[i]]$focal$multilocus)] <- amDatasetFocal$missingCode

      pairwiseMatches[[i]]$match <- list(
        index = NULL,
        metaData = NULL,
        multilocus = NULL,
        score = NULL
      )

      if (length(pairwiseMatchesWhich) == 0) {
        pairwiseMatches[[i]]$match$index <- "None"
        pairwiseMatches[[i]]$match$metaData <- if (!is.null(metaDataFocal[i]))
          ""
        else
          NULL
        pairwiseMatches[[i]]$match$multilocus <- matrix("", 1, length(columnNames))
        pairwiseMatches[[i]]$match$score <- ""
        pairwiseMatches[[i]]$match$flags <- matrix(1, nrow = 1, ncol = length(columnNames))
        pairwiseMatches[[i]]$match$perfect <- 0
        pairwiseMatches[[i]]$match$partial <- 0
      } else {
        pairwiseMatches[[i]]$match$index <- indexComparison[pairwiseMatchesWhich]
        pairwiseMatches[[i]]$match$metaData <- metaDataComparison[pairwiseMatchesWhich]
        pairwiseMatches[[i]]$match$multilocus <- comparisonGenotypes[pairwiseMatchesWhich, , drop = FALSE]
        pairwiseMatches[[i]]$match$score <- format(signif(pairwiseMatchesScores, 2))

        if (nrow(pairwiseMatches[[i]]$match$multilocus) > 1) {
          ord <- order(pairwiseMatchesScores, decreasing = TRUE)
          pairwiseMatches[[i]]$match$index <- pairwiseMatches[[i]]$match$index[ord]
          pairwiseMatches[[i]]$match$metaData <- pairwiseMatches[[i]]$match$metaData[ord]
          pairwiseMatches[[i]]$match$multilocus <- pairwiseMatches[[i]]$match$multilocus[ord, , drop = FALSE]
          pairwiseMatches[[i]]$match$score <- pairwiseMatches[[i]]$match$score[ord]
        }

        pairwiseMatches[[i]]$match$flags <- matrix(
          1,
          nrow(pairwiseMatches[[i]]$match$multilocus),
          ncol(pairwiseMatches[[i]]$match$multilocus)
        )
        focal_vec <- as.vector(pairwiseMatches[[i]]$focal$multilocus)
        mismatch_mask <- comparisonGenotypes[pairwiseMatchesWhich, , drop = FALSE] != focal_vec[col(comparisonGenotypes[pairwiseMatchesWhich, , drop = FALSE])]
        pairwiseMatches[[i]]$match$flags[mismatch_mask] <- 0
        pairwiseMatches[[i]]$match$flags[is.na(pairwiseMatches[[i]]$match$multilocus)] <- 2
        pairwiseMatches[[i]]$match$flags[, which(pairwiseMatches[[i]]$focal$flags == 2)] <- 2

        pairwiseMatches[[i]]$match$perfect <- sum(pairwiseMatchesScores == 1)
        pairwiseMatches[[i]]$match$partial <- sum(pairwiseMatchesScores < 1)
      }

      dimnames(pairwiseMatches[[i]]$match$multilocus)[[2]] <- columnNames
      pairwiseMatches[[i]]$match$multilocus[is.na(pairwiseMatches[[i]]$match$multilocus)] <- amDatasetFocal$missingCode
    }

    amPairwise <- list()
    amPairwise$pairwise <- pairwiseMatches
    amPairwise$missingCode <- amDatasetFocal$missingCode
    amPairwise$matchThreshold <- matchThreshold
    amPairwise$alleleMismatch <- alleleMismatch
    amPairwise$missingMethod <- missingMethod
    amPairwise$focalDatasetN <- nrow(amDatasetFocal$multilocus)
    amPairwise$comparisonDatasetN <- nrow(amDatasetComparison$multilocus)
    amPairwise$focalIsComparison <- identical(dim(amDatasetFocal$multilocus),
                                              dim(amDatasetComparison$multilocus)) &&
      all(amDatasetFocal$multilocus == amDatasetComparison$multilocus)

    class(amPairwise) <- "amPairwise"
    return(amPairwise)
  }


#### summary.amPairwise() ####
summary.amPairwise <- function(object,
                               html = NULL,
                               csv = NULL,
                               ...) {
  if (!inherits(object, "amPairwise"))
    stop("allelematch: this function requires an \"amPairwise\" object",
         call. = FALSE)

  if (!is.null(html)) {
    if (is.logical(html) &&
        (html == TRUE))
      amHTML.amPairwise(object)
    else if (is.logical(html) &&
             (html == FALSE))
      html <- NULL
    else
      amHTML.amPairwise(object, htmlFile = html)
  }
  if (!is.null(csv))
    amCSV.amPairwise(object, csvFile = csv)

  if (is.null(html) && is.null(csv)) {
    comparison_str <- if (object$focalIsComparison)
      "focal dataset compared against itself\n"
    else
      paste0("comparison dataset N=", object$comparisonDatasetN, "\n")

    cat(
      "allelematch\npairwise analysis\n\n",
      "focal dataset N=",
      object$focalDatasetN,
      "\n",
      comparison_str,
      "missing data represented by: ",
      object$missingCode,
      "\n",
      "missing data matching method: ",
      object$missingMethod,
      "\n",
      "alleleMismatch (m-hat; maximum number of mismatching alleles): ",
      object$alleleMismatch,
      "\n",
      "matchThreshold (s-hat; lowest matching score returned): ",
      object$matchThreshold,
      "\n\n",
      "score flags:\n",
      "*101 allele does not match\n",
      "+101 allele is missing\n\n",
      sep = ""
    )

    y <- object$pairwise
    for (i in 1:length(y)) {
      cat("(", i, " of ", length(y), ")\n", sep = "")
      if (is.null(y[[i]]$focal$metaData))
        y[[i]]$focal$metaData <- ""
      if (is.null(y[[i]]$match$metaData))
        y[[i]]$match$metaData <- rep("", length(y[[i]]$match$index))

      showFocalFlags <- y[[i]]$focal$multilocus
      showFocalFlags[y[[i]]$focal$flags == 2] <- "+"
      showFocalFlags[y[[i]]$focal$flags == 0] <- "*"
      showFocalFlags[y[[i]]$focal$flags == 1] <- ""
      showFocal <- matrix(
        paste0(showFocalFlags, y[[i]]$focal$multilocus),
        nrow(showFocalFlags),
        ncol(showFocalFlags)
      )
      dimnames(showFocal) <- dimnames(showFocalFlags)

      showMatchFlags <- y[[i]]$match$multilocus
      showMatchFlags[y[[i]]$match$flags == 2] <- "+"
      showMatchFlags[y[[i]]$match$flags == 0] <- "*"
      showMatchFlags[y[[i]]$match$flags == 1] <- ""
      showMatch <- matrix(
        paste0(showMatchFlags, y[[i]]$match$multilocus),
        nrow(showMatchFlags),
        ncol(showMatchFlags)
      )
      dimnames(showMatch) <- dimnames(showMatchFlags)

      print(data.frame(
        rbind(
          data.frame(showFocal, score = ""),
          data.frame(showMatch, score = y[[i]]$match$score)
        ),
        row.names = paste0(c("FOCAL ", rep(
          "MATCH ", length(y[[i]]$match$index)
        )), format(c(
          y[[i]]$focal$index, y[[i]]$match$index
        )), " ", format(
          c(y[[i]]$focal$metaData, y[[i]]$match$metaData)
        ))
      ))
      cat(
        y[[i]]$match$perfect,
        " perfect matches found. ",
        y[[i]]$match$partial,
        " partial matches found.\n\n\n",
        sep = ""
      )
    }
  }
}


#### amCluster() ####
amCluster <-
  function(amDatasetFocal,
           runUntilSingletons = TRUE,
           cutHeight = 0.3,
           missingMethod = 2,
           consensusMethod = 1,
           clusterMethod = "complete") {
    if (!(class(amDatasetFocal) %in% c("amDataset", "amInterpolate", "amCluster"))) {
      stop(
        "allelematch: amDatasetFocal must be an object of class \"amDataset\" or for recursive use, class \"amCluster\"",
        call. = FALSE
      )
    }
    if (inherits(amDatasetFocal, "amCluster"))
      amDatasetFocal <- amDatasetFocal$unique

    if (inherits(amDatasetFocal, "amInterpolate")) {
      reClass <- amDatasetFocal
      amDatasetFocal <- list(
        index = reClass$index,
        metaData = reClass$metaData,
        multilocus = reClass$multilocus,
        missingCode = reClass$missingCode
      )
      class(amDatasetFocal) <- "amDataset"
    }

    originalFocalDatasetN <- nrow(amDatasetFocal$multilocus)
    if (!(missingMethod %in% c(1, 2)))
      stop("allelematch: missingMethod must equal 1 or 2", call. = FALSE)
    if (!(consensusMethod %in% c(1, 2, 3, 4)))
      stop("allelematch: consensusMethod must equal 1, 2, 3 or 4", call. = FALSE)
    if (!(tolower(clusterMethod) %in% c("complete", "average", "single")))
      stop("allelematch: clusterMethod must be \"complete\" or \"average\"",
           call. = FALSE)

    totalRuns <- 0
    repeat {
      totalRuns <- totalRuns + 1
      if (is.null(dim(amDatasetFocal$multilocus)))
        break

      dissimMatrix <- amMatrix(amDatasetFocal, missingMethod = missingMethod)

      tryCatch(
        dendro <- stats::hclust(stats::as.dist(dissimMatrix), method = clusterMethod),
        error = function(x)
          stop(
            "allelematch: error in clustering, try a different clusterMethod",
            call. = FALSE
          )
      )

      tryCatch(
        uniqueIndex <- dynamicTreeCut::cutreeHybrid(
          dendro,
          dissimMatrix,
          cutHeight = cutHeight,
          minClusterSize = 1,
          verbose = 0
        )$labels + 1,
        error = function(x)
          stop(
            "allelematch: error in dynamic tree cutting; several causes for this error; have you installed dynamicTreeCut package?  install.packages(\"dynamicTreeCut\")",
            call. = FALSE
          )
      )
      numLabels <- length(unique(uniqueIndex))

      clusterAnalysis <- list()
      clusterAnalysis$cluster <- vector("list", numLabels - 1)

      j <- 0
      for (i in 1:numLabels) {
        thisGenotype <- amDatasetFocal$multilocus[uniqueIndex == i, , drop = FALSE]
        dimnames(thisGenotype) <- list(NULL, dimnames(amDatasetFocal$multilocus)[[2]])
        thisIndex <- amDatasetFocal$index[uniqueIndex == i]
        thisMetaData <- amDatasetFocal$metaData[uniqueIndex == i]

        if ((i == 1) && (length(thisGenotype) > 1)) {
          modifiedDissimMatrix <- dissimMatrix
          diag(modifiedDissimMatrix) <- 1

          maxScoreSingletons <- apply(matrix(
            1 - modifiedDissimMatrix[uniqueIndex == 1, ],
            sum(uniqueIndex == 1),
            ncol(modifiedDissimMatrix)
          ), 1, max)
          maxScoreSingletonsWhich <- apply(matrix(
            1 - modifiedDissimMatrix[uniqueIndex == 1, ],
            sum(uniqueIndex == 1),
            ncol(modifiedDissimMatrix)
          ),
          1,
          which.max)
          clusterAnalysis$singletons <- vector("list", sum(uniqueIndex == 1))

          if (sum(uniqueIndex == 1) > 1) {
            reorderSingletons <- order(maxScoreSingletons, decreasing = TRUE)
            maxScoreSingletons <- maxScoreSingletons[reorderSingletons]
            maxScoreSingletonsWhich <- maxScoreSingletonsWhich[reorderSingletons]
            thisGenotype <- thisGenotype[reorderSingletons, ]
            thisIndex <- thisIndex[reorderSingletons]
            thisMetaData <- thisMetaData[reorderSingletons]
          }

          for (k in 1:sum(uniqueIndex == 1)) {
            clusterAnalysis$singletons[[k]]$focal$index <- thisIndex[k]
            clusterAnalysis$singletons[[k]]$focal$metaData <- thisMetaData[k]
            clusterAnalysis$singletons[[k]]$focal$multilocus <- t(thisGenotype[k, ])
            clusterAnalysis$singletons[[k]]$focal$flags <- matrix(
              1,
              nrow(clusterAnalysis$singletons[[k]]$focal$multilocus),
              ncol(clusterAnalysis$singletons[[k]]$focal$multilocus)
            )
            clusterAnalysis$singletons[[k]]$focal$flags[clusterAnalysis$singletons[[k]]$focal$multilocus == amDatasetFocal$missingCode] <- 2

            clusterAnalysis$singletons[[k]]$match$index <- amDatasetFocal$index[maxScoreSingletonsWhich[k]]
            clusterAnalysis$singletons[[k]]$match$metaData <- amDatasetFocal$metaData[maxScoreSingletonsWhich[k]]
            clusterAnalysis$singletons[[k]]$match$multilocus <- t(amDatasetFocal$multilocus[maxScoreSingletonsWhich[k], ])
            clusterAnalysis$singletons[[k]]$match$score <- as.character(signif(maxScoreSingletons[k], 2))
            clusterAnalysis$singletons[[k]]$match$flags <- matrix(
              1,
              nrow(clusterAnalysis$singletons[[k]]$match$multilocus),
              ncol(clusterAnalysis$singletons[[k]]$match$multilocus)
            )

            clusterAnalysis$singletons[[k]]$match$flags[matrix(
              clusterAnalysis$singletons[[k]]$focal$multilocus,
              nrow(clusterAnalysis$singletons[[k]]$match$multilocus),
              ncol(clusterAnalysis$singletons[[k]]$match$multilocus),
              byrow = TRUE
            ) != clusterAnalysis$singletons[[k]]$match$multilocus] <- 0
            clusterAnalysis$singletons[[k]]$match$flags[clusterAnalysis$singletons[[k]]$match$multilocus == amDatasetFocal$missingCode] <- 2
            clusterAnalysis$singletons[[k]]$match$flags[, which(clusterAnalysis$singletons[[k]]$focal$flags == 2)] <- 2
          }
        } else if (i == 1) {
          clusterAnalysis$singletons <- list()
        } else {
          j <- j + 1
          score <- amMatrix(
            amDataset(
              cbind(1:nrow(thisGenotype), thisGenotype),
              indexColumn = 1,
              missingCode = amDatasetFocal$missingCode
            ),
            missingMethod = missingMethod
          )

          if (consensusMethod %in% c(1, 2)) {
            lowerTriScore <- score
            lowerTriScore[upper.tri(score, diag = TRUE)] <- 0
            consensusIndex <- which.min(rowSums(lowerTriScore))
            clusterAnalysis$cluster[[j]]$focal$index <- thisIndex[consensusIndex]
            clusterAnalysis$cluster[[j]]$focal$metaData <- thisMetaData[consensusIndex]
            clusterAnalysis$cluster[[j]]$focal$multilocus <- t(thisGenotype[consensusIndex, ])

            if (consensusMethod == 2 &&
                sum(
                  clusterAnalysis$cluster[[j]]$focal$multilocus == amDatasetFocal$missingCode
                ) > 0) {
              modes <- apply(thisGenotype, 2, function(x) {
                modeAllele <- attr(sort(table(x), decreasing = TRUE), "name")
                if (length(modeAllele) > 1 &&
                    modeAllele[1] == amDatasetFocal$missingCode)
                  modeAllele[2]
                else
                  modeAllele[1]
              })
              reassignThese <- clusterAnalysis$cluster[[j]]$focal$multilocus == amDatasetFocal$missingCode
              clusterAnalysis$cluster[[j]]$focal$multilocus[reassignThese] <- modes[reassignThese]
              interpolated <- reassignThese &
                (modes != amDatasetFocal$missingCode)
            }
          } else if (consensusMethod %in% c(3, 4)) {
            missingAlleles <- rowSums(thisGenotype == amDatasetFocal$missingCode)
            consensusIndex <- which.min(missingAlleles)
            clusterAnalysis$cluster[[j]]$focal$index <- thisIndex[consensusIndex]
            clusterAnalysis$cluster[[j]]$focal$metaData <- thisMetaData[consensusIndex]
            clusterAnalysis$cluster[[j]]$focal$multilocus <- t(thisGenotype[consensusIndex, ])

            if (consensusMethod == 4 &&
                sum(
                  clusterAnalysis$cluster[[j]]$focal$multilocus == amDatasetFocal$missingCode
                ) > 0) {
              modes <- apply(thisGenotype, 2, function(x) {
                modeAllele <- attr(sort(table(x), decreasing = TRUE), "name")
                if (length(modeAllele) > 1 &&
                    modeAllele[1] == amDatasetFocal$missingCode)
                  modeAllele[2]
                else
                  modeAllele[1]
              })
              reassignThese <- clusterAnalysis$cluster[[j]]$focal$multilocus == amDatasetFocal$missingCode
              clusterAnalysis$cluster[[j]]$focal$multilocus[reassignThese] <- modes[reassignThese]
              interpolated <- reassignThese &
                (modes != amDatasetFocal$missingCode)
            }
          }

          clusterAnalysis$cluster[[j]]$focal$flags <- matrix(
            1,
            nrow(clusterAnalysis$cluster[[j]]$focal$multilocus),
            ncol(clusterAnalysis$cluster[[j]]$focal$multilocus)
          )
          clusterAnalysis$cluster[[j]]$focal$flags[clusterAnalysis$cluster[[j]]$focal$multilocus == amDatasetFocal$missingCode] <- 2
          if (exists("interpolated"))
            clusterAnalysis$cluster[[j]]$focal$flags[interpolated] <- 3

          reorderMatch <- order(1 - score[consensusIndex, ], decreasing = TRUE)
          clusterAnalysis$cluster[[j]]$match$index <- thisIndex[reorderMatch]
          clusterAnalysis$cluster[[j]]$match$metaData <- thisMetaData[reorderMatch]
          clusterAnalysis$cluster[[j]]$match$multilocus <- thisGenotype[reorderMatch, ]
          clusterAnalysis$cluster[[j]]$match$score <- as.character(signif(1 - score[consensusIndex, ], 2)[reorderMatch])

          clusterAnalysis$cluster[[j]]$match$flags <- matrix(
            1,
            nrow(clusterAnalysis$cluster[[j]]$match$multilocus),
            ncol(clusterAnalysis$cluster[[j]]$match$multilocus)
          )
          clusterAnalysis$cluster[[j]]$match$flags[matrix(
            clusterAnalysis$cluster[[j]]$focal$multilocus,
            nrow(clusterAnalysis$cluster[[j]]$match$multilocus),
            ncol(clusterAnalysis$cluster[[j]]$match$multilocus),
            byrow = TRUE
          ) != clusterAnalysis$cluster[[j]]$match$multilocus] <- 0
          clusterAnalysis$cluster[[j]]$match$flags[clusterAnalysis$cluster[[j]]$match$multilocus == amDatasetFocal$missingCode] <- 2
          clusterAnalysis$cluster[[j]]$match$flags[, which(clusterAnalysis$cluster[[j]]$focal$flags == 2)] <- 2

          clusterAnalysis$cluster[[j]]$match$perfect <- sum(apply(clusterAnalysis$cluster[[j]]$match$flags, 1, function(x)
            all(x == 1)))
          clusterAnalysis$cluster[[j]]$match$partial <- nrow(clusterAnalysis$cluster[[j]]$match$flags) - clusterAnalysis$cluster[[j]]$match$perfect
        }
      }

      if (length(clusterAnalysis$cluster) == 0) {
        clusterAnalysis$unique$index <- do.call(c,
                                                lapply(clusterAnalysis$singletons, function(x)
                                                  x$focal$index))
        clusterAnalysis$unique$metaData <- do.call(c,
                                                   lapply(clusterAnalysis$singletons, function(x)
                                                     x$focal$metaData))
        clusterAnalysis$unique$multilocus <- do.call(rbind,
                                                     lapply(clusterAnalysis$singletons, function(x)
                                                       x$focal$multilocus))
        clusterAnalysis$unique$uniqueType <- rep("SINGLETON", length(clusterAnalysis$singletons))
      } else if (length(clusterAnalysis$singletons) == 0) {
        clusterAnalysis$unique$index <- do.call(c,
                                                lapply(clusterAnalysis$cluster, function(x)
                                                  x$focal$index))
        clusterAnalysis$unique$metaData <- do.call(c,
                                                   lapply(clusterAnalysis$cluster, function(x)
                                                     x$focal$metaData))
        clusterAnalysis$unique$multilocus <- do.call(rbind,
                                                     lapply(clusterAnalysis$cluster, function(x)
                                                       x$focal$multilocus))
        clusterAnalysis$unique$uniqueType <- rep("CONSENSUS", length(clusterAnalysis$cluster))
      } else {
        clusterAnalysis$unique$index <- c(do.call(
          c,
          lapply(clusterAnalysis$cluster, function(x)
            x$focal$index)
        ), do.call(
          c,
          lapply(clusterAnalysis$singletons, function(x)
            x$focal$index)
        ))
        clusterAnalysis$unique$metaData <- c(do.call(
          c,
          lapply(clusterAnalysis$cluster, function(x)
            x$focal$metaData)
        ), do.call(
          c,
          lapply(clusterAnalysis$singletons, function(x)
            x$focal$metaData)
        ))
        clusterAnalysis$unique$multilocus <- rbind(do.call(
          rbind,
          lapply(clusterAnalysis$cluster, function(x)
            x$focal$multilocus)
        ), do.call(
          rbind,
          lapply(clusterAnalysis$singletons, function(x)
            x$focal$multilocus)
        ))
        clusterAnalysis$unique$uniqueType <- c(rep("CONSENSUS", length(clusterAnalysis$cluster)), rep("SINGLETON", length(clusterAnalysis$singletons)))
      }

      if (sum(is.na(suppressWarnings(
        as.numeric(clusterAnalysis$unique$index)
      ))) > 0) {
        orderUnique <- order(clusterAnalysis$unique$index)
      } else {
        orderUnique <- order(as.numeric(clusterAnalysis$unique$index))
      }
      clusterAnalysis$unique$index <- clusterAnalysis$unique$index[orderUnique]
      clusterAnalysis$unique$metaData <- clusterAnalysis$unique$metaData[orderUnique]
      clusterAnalysis$unique$multilocus <- clusterAnalysis$unique$multilocus[orderUnique, ]
      clusterAnalysis$unique$uniqueType <- clusterAnalysis$unique$uniqueType[orderUnique]
      clusterAnalysis$unique$missingCode <- amDatasetFocal$missingCode
      class(clusterAnalysis$unique) <- "amDataset"

      clusterAnalysis$cutHeight <- cutHeight
      clusterAnalysis$consensusMethod <- consensusMethod
      clusterAnalysis$missingMethod <- missingMethod
      clusterAnalysis$clusterMethod <- clusterMethod
      clusterAnalysis$missingCode <- amDatasetFocal$missingCode
      clusterAnalysis$focalDatasetN <- originalFocalDatasetN
      clusterAnalysis$totalRuns <- totalRuns
      clusterAnalysis$runUntilSingletons <- runUntilSingletons
      class(clusterAnalysis) <- "amCluster"

      if (runUntilSingletons &&
          length(clusterAnalysis$cluster) > 0)
        amDatasetFocal <- clusterAnalysis$unique
      else
        break
    }

    if (is.null(clusterAnalysis$unique$multilocus))
      stop(
        "allelematch: amCluster: no clusters formed. Please set cutHeight lower and run again",
        call. = FALSE
      )
    if (is.null(dim(clusterAnalysis$unique$multilocus))) {
      tmpMultilocus <- matrix(
        clusterAnalysis$unique$multilocus,
        1,
        length(clusterAnalysis$unique$multilocus),
        byrow = FALSE
      )
      dimnames(tmpMultilocus)[[2]] <- names(clusterAnalysis$unique$multilocus)
      clusterAnalysis$unique$multilocus <- tmpMultilocus
    }
    return(clusterAnalysis)
  }


#### summary.amCluster() ####
summary.amCluster <- function(object,
                              html = NULL,
                              csv = NULL,
                              ...) {
  if (!inherits(object, "amCluster"))
    stop("allelematch: this function requires an \"amCluster\" object",
         call. = FALSE)

  if (!is.null(html)) {
    if (is.logical(html) &&
        (html == TRUE))
      amHTML.amCluster(object)
    else if (is.logical(html) &&
             (html == FALSE))
      html <- NULL
    else
      amHTML.amCluster(object, htmlFile = html)
  }
  if (!is.null(csv))
    amCSV.amCluster(object, csvFile = csv)

  if (is.null(html) && is.null(csv)) {
    cat(
      "allelematch\namCluster object\n\n",
      "Focal dataset N=",
      object$focalDatasetN,
      "\n",
      "Unique genotypes (total): ",
      nrow(object$unique$multilocus),
      "\n",
      "Unique genotypes (by cluster consensus): ",
      length(object$cluster),
      "\n",
      "Unique genotypes (singletons): ",
      length(object$singletons),
      "\n\n",
      "Missing data represented by: ",
      object$missingCode,
      "\n",
      "Missing data matching method: ",
      object$missingMethod,
      "\n",
      "Clustered genotypes consensus method: ",
      object$consensusMethod,
      "\n",
      "Hierarchical clustering method: ",
      object$clusterMethod,
      "\n",
      "Dynamic tree cutting height (cutHeight): ",
      object$cutHeight,
      "\n\n",
      "Run until only singletons: ",
      object$runUntilSingletons,
      "\n",
      "Runs: ",
      object$totalRuns,
      "\n\n",
      "Score flags:\n",
      "*101 Allele does not match\n",
      "+101 Allele is missing\n",
      "[101] Allele is interpolated in consensus from non-missing cluster members (consensusMethod = 2, 4 only)\n\n",
      sep = ""
    )

    y <- object$cluster
    if (length(y) > 0) {
      for (i in 1:length(y)) {
        cat("(Cluster ", i, " of ", length(y), ")\n", sep = "")
        if (is.null(y[[i]]$focal$metaData))
          y[[i]]$focal$metaData <- ""
        if (is.null(y[[i]]$match$metaData))
          y[[i]]$match$metaData <- rep("", length(y[[i]]$match$index))

        showFocalFlags <- y[[i]]$focal$multilocus
        showFocalFlags[y[[i]]$focal$flags == 2] <- "+"
        showFocalFlags[y[[i]]$focal$flags == 0] <- "*"
        showFocalFlags[y[[i]]$focal$flags == 1] <- ""
        showFocalFlags[y[[i]]$focal$flags == 3] <- "["
        showFocal <- matrix(
          paste0(showFocalFlags, y[[i]]$focal$multilocus),
          nrow(showFocalFlags),
          ncol(showFocalFlags)
        )
        showFocal[, grepl("\\[", showFocal)] <- paste0(showFocal[, grepl("\\[", showFocal)], "]")
        dimnames(showFocal) <- dimnames(showFocalFlags)

        showMatchFlags <- y[[i]]$match$multilocus
        showMatchFlags[y[[i]]$match$flags == 2] <- "+"
        showMatchFlags[y[[i]]$match$flags == 0] <- "*"
        showMatchFlags[y[[i]]$match$flags == 1] <- ""
        showMatch <- matrix(
          paste0(showMatchFlags, y[[i]]$match$multilocus),
          nrow(showMatchFlags),
          ncol(showMatchFlags)
        )
        dimnames(showMatch) <- dimnames(showMatchFlags)

        print(data.frame(
          rbind(
            data.frame(showFocal, score = ""),
            data.frame(showMatch, score = y[[i]]$match$score)
          ),
          row.names = paste0(c(
            "CONSENSUS ", rep("MEMBER    ", length(y[[i]]$match$index))
          ), format(
            c(y[[i]]$focal$index, y[[i]]$match$index)
          ), " ", format(
            c(y[[i]]$focal$metaData, y[[i]]$match$metaData)
          ))
        ))
        cat(
          y[[i]]$match$perfect,
          " perfect matches found. ",
          y[[i]]$match$partial,
          " partial matches found.\n\n\n",
          sep = ""
        )
      }
      cat("(Clusters END)\n\n\n")
    }

    y <- object$singletons
    if (length(y) > 0) {
      for (i in 1:length(y)) {
        cat("(Singleton ", i, " of ", length(y), ")\n", sep = "")
        if (is.null(y[[i]]$focal$metaData))
          y[[i]]$focal$metaData <- ""
        if (is.null(y[[i]]$match$metaData))
          y[[i]]$match$metaData <- rep("", length(y[[i]]$match$index))

        showFocalFlags <- y[[i]]$focal$multilocus
        showFocalFlags[y[[i]]$focal$flags == 2] <- "+"
        showFocalFlags[y[[i]]$focal$flags == 0] <- "*"
        showFocalFlags[y[[i]]$focal$flags == 1] <- ""
        showFocal <- matrix(
          paste0(showFocalFlags, y[[i]]$focal$multilocus),
          nrow(showFocalFlags),
          ncol(showFocalFlags)
        )
        dimnames(showFocal) <- dimnames(showFocalFlags)

        showMatchFlags <- y[[i]]$match$multilocus
        showMatchFlags[y[[i]]$match$flags == 2] <- "+"
        showMatchFlags[y[[i]]$match$flags == 0] <- "*"
        showMatchFlags[y[[i]]$match$flags == 1] <- ""
        showMatch <- matrix(
          paste0(showMatchFlags, y[[i]]$match$multilocus),
          nrow(showMatchFlags),
          ncol(showMatchFlags)
        )
        dimnames(showMatch) <- dimnames(showMatchFlags)

        print(data.frame(
          rbind(
            data.frame(showFocal, score = ""),
            data.frame(showMatch, score = y[[i]]$match$score)
          ),
          row.names = paste0(c(
            "SINGLETON ", rep("CLOSEST   ", length(y[[i]]$match$index))
          ), format(
            c(y[[i]]$focal$index, y[[i]]$match$index)
          ), " ", format(
            c(y[[i]]$focal$metaData, y[[i]]$match$metaData)
          ))
        ))
        cat("\n\n")
      }
      cat("(Singletons END)\n\n\n")
    }

    y <- object$unique
    cat("(Unique genotypes)\n")
    if (is.null(y$metaData))
      y$metaData <- ""
    print(data.frame(
      y$multilocus,
      row.names = paste0(y$uniqueType, " ", format(y$index), " ", format(y$metaData))
    ))
  }
}


#### amAlleleFreq() ####
amAlleleFreq <- function(amDatasetFocal, multilocusMap = NULL) {
  if (!inherits(amDatasetFocal, "amDataset"))
    stop("allelematch: amDatasetFocal must be an object of class \"amDataset\"",
         call. = FALSE)

  if (is.null(multilocusMap)) {
    if ((ncol(amDatasetFocal$multilocus) %% 2) != 0) {
      stop(
        "allelematch: there are an odd number of genotype columns in amDatasetFocal; Please specify multilocusMap manually",
        call. = FALSE
      )
    } else {
      cat(
        "allelematch: assuming genotype columns are in pairs, representing ",
        ncol(amDatasetFocal$multilocus) / 2,
        " loci\n",
        sep = ""
      )
    }
    multilocusMap <- rep(1:(ncol(amDatasetFocal$multilocus) / 2), each = 2)
  } else if (length(multilocusMap) != ncol(amDatasetFocal$multilocus)) {
    stop(
      "allelematch: multilocusMap must be a vector of integers or strings giving the mappings onto loci for all genotype columns in amDatasetFocal;\n",
      "             Example: sex followed by 4 diploid loci in paired columns could be coded: mutlilocusMap=c(1,2,2,3,3,4,4,5,5)\n",
      "             or as: multilocusMap=c(\"SEX\",\"LOC1\",\"LOC1\",\"LOC2\",\"LOC2\",\"LOC3\",\"LOC3\",\"LOC4\",\"LOC4\")",
      call. = FALSE
    )
  } else if (sum(table(multilocusMap) > 2) > 0) {
    stop(
      "allelematch: multilocusMap indicates that a locus occurs in three or more columns;  this situation is not yet handled",
      call. = FALSE
    )
  }
  multilocusMap <- as.integer(as.factor(multilocusMap))

  alleleFreq <- list(multilocusMap = multilocusMap, loci = vector("list", max(multilocusMap)))

  for (locus in 1:max(multilocusMap)) {
    thisLocusIndex <- which(multilocusMap == locus)
    alleleFreq$loci[[locus]]$name <- paste(dimnames(amDatasetFocal$multilocus)[[2]][thisLocusIndex], collapse = "-")
    alleleFreq$loci[[locus]]$columnNames <- dimnames(amDatasetFocal$multilocus)[[2]][thisLocusIndex]
    thisLocus <- as.character(amDatasetFocal$multilocus[, thisLocusIndex])
    thisLocus[thisLocus == amDatasetFocal$missingCode] <- NA
    thisLocusUnique <- unique(thisLocus)[!is.na(unique(thisLocus))]
    alleleFreq$loci[[locus]]$alleleFreq <- sort(sapply(thisLocusUnique, function(x)
      sum(thisLocus == x, na.rm = TRUE) / length(thisLocus[!is.na(thisLocus)])),
      decreasing = TRUE)
    alleleFreq$loci[[locus]]$missingFreq <- sum(is.na(thisLocus)) / length(thisLocus)
    alleleFreq$loci[[locus]]$numAlleles <- length(thisLocusUnique)
    # Fix: Calculate PIC only if there are at least 2 alleles
    if (alleleFreq$loci[[locus]]$numAlleles >= 2) {
      alleleFreq$loci[[locus]]$PIC <- 1 - sum(alleleFreq$loci[[locus]]$alleleFreq^2) -
        sum(apply(t(combn(
          1:length(alleleFreq$loci[[locus]]$alleleFreq), 2
        )), 1, function(x)
          prod(alleleFreq$loci[[locus]]$alleleFreq[x]^2)))
    } else {
      # For monomorphic loci, PIC is defined as 0
      alleleFreq$loci[[locus]]$PIC <- 0
    }
  }
  class(alleleFreq) <- "amAlleleFreq"
  return(alleleFreq)
}


#### print.amAlleleFreq() ####
print.amAlleleFreq <- function(x, ...) {
  cat(
    "allelematch\namAlleleFreq object\nFrequencies calculated after removal of missing data\n"
  )
  y <- x$loci
  for (i in 1:length(y)) {
    cat("\n", y[[i]]$name, " (", y[[i]]$numAlleles, " alleles)\n", sep = "")
    for (j in 1:length(y[[i]]$alleleFreq)) {
      cat(
        "\tAllele\t",
        names(y[[i]]$alleleFreq)[j],
        "\t",
        signif(y[[i]]$alleleFreq[j], 3),
        "\n",
        sep = ""
      )
    }
  }
}


#### amUnique() ####
amUnique <-
  function(amDatasetFocal,
           multilocusMap = NULL,
           alleleMismatch = NULL,
           matchThreshold = NULL,
           cutHeight = NULL,
           doPsib = "missing",
           consensusMethod = 1,
           verbose = TRUE) {
    if (!inherits(amDatasetFocal, "amDataset"))
      stop("allelematch: amDatasetFocal must be an object of class \"amDataset\"",
           call. = FALSE)

    if (is.null(multilocusMap)) {
      if ((ncol(amDatasetFocal$multilocus) %% 2) != 0) {
        stop(
          "allelematch: there are an odd number of genotype columns in amDatasetFocal; Please specify multilocusMap manually",
          call. = FALSE
        )
      } else if (verbose) {
        cat(
          "allelematch: assuming genotype columns are in pairs, representing ",
          ncol(amDatasetFocal$multilocus) / 2,
          " loci\n",
          sep = ""
        )
      }
      multilocusMap <- rep(1:(ncol(amDatasetFocal$multilocus) / 2), each = 2)
    } else if (length(multilocusMap) != ncol(amDatasetFocal$multilocus)) {
      stop(
        "allelematch: multilocusMap must be a vector of integers or strings giving the mappings onto loci for all genotype columns in amDatasetFocal;
             Example: sex followed by 4 diploid loci in paired columns could be coded: mutlilocusMap=c(1,2,2,3,3,4,4,5,5)
             or as: multilocusMap=c(\"GENDER\",\"LOC1\",\"LOC1\",\"LOC2\",\"LOC2\",\"LOC3\",\"LOC3\",\"LOC4\",\"LOC4\")",
        call. = FALSE
      )
    } else if (sum(table(multilocusMap) > 2) > 0) {
      stop(
        "allelematch: multilocusMap indicates that a locus occurs in three or more columns; this situation is not yet handled",
        call. = FALSE
      )
    }
    multilocusMap <- as.integer(as.factor(multilocusMap))

    if (sum(!(c(
      is.null(alleleMismatch),
      is.null(matchThreshold),
      is.null(cutHeight)
    ))) != 1)
      stop("allelematch: please specify alleleMismatch OR matchThreshold OR cutHeight",
           call. = FALSE)
    if (length(c(alleleMismatch, matchThreshold, cutHeight)) > 1)
      stop(
        "allelematch: please provide a single parameter value for alleleMismatch OR matchThreshold OR cutHeight. Use amUniqueProfile() to examine a range of values",
        call. = FALSE
      )

    if (!is.null(matchThreshold)) {
      if ((matchThreshold < 0) ||
          (matchThreshold > 1))
        stop("allelematch: matchThreshold must be between 0 and 1", call. = FALSE)
      cutHeight <- 1 - matchThreshold
      alleleMismatch <- round((1 - matchThreshold) * length(multilocusMap), 2)
    } else if (!is.null(alleleMismatch)) {
      matchThreshold <- 1 - (alleleMismatch / length(multilocusMap))
      cutHeight <- 1 - matchThreshold
    } else if (!is.null(cutHeight)) {
      if ((cutHeight < 0) ||
          (cutHeight > 1))
        stop("allelematch: cutHeight must be greater than 0 and less than 1",
             call. = FALSE)
      matchThreshold <- 1 - cutHeight
      alleleMismatch <- round((1 - matchThreshold) * length(multilocusMap), 2)
    }

    if (matchThreshold == 1 && cutHeight == 0) {
      if (verbose)
        cat(
          "allelematch: cutHeight cannot be zero. Setting cutHeight=0.00001. This will return perfect matches\n"
        )
      cutHeight <- 0.00001
    }

    if (verbose)
      cat("allelematch: amUnique: Clustering genotypes\n")
    clusterAnalysis <- amCluster(
      amDatasetFocal,
      cutHeight = cutHeight,
      runUntilSingletons = TRUE,
      consensusMethod = consensusMethod
    )

    if (verbose)
      cat(
        "allelematch: amUnique: Comparing unique genotypes identified by clustering to all samples\n"
      )
    clusterAnalysisPairwise <- amPairwise(clusterAnalysis$unique, amDatasetFocal, matchThreshold = matchThreshold)

    if (verbose)
      cat(
        "allelematch: amUnique: Determining allele frequencies of unique genotypes identified by cluster\n"
      )
    clusterAnalysisAlleleFreq <- amAlleleFreq(clusterAnalysis$unique, multilocusMap = multilocusMap)

    if (verbose)
      cat("allelematch: amUnique: Finding Psib\n")
    uniqueAnalysis <- clusterAnalysisPairwise
    class(uniqueAnalysis) <- "amUnique"
    for (i in 1:length(uniqueAnalysis$pairwise)) {
      uniqueAnalysis$pairwise[[i]]$match$psib <- rep(NA, nrow(uniqueAnalysis$pairwise[[i]]$match$multilocus))
      uniqueAnalysis$pairwise[[i]]$match$rowFlag <- rep("", nrow(uniqueAnalysis$pairwise[[i]]$match$multilocus))
    }

    multilocusMap <- clusterAnalysisAlleleFreq$multilocusMap
    for (iPairwise in 1:length(uniqueAnalysis$pairwise)) {
      thisPairwise <- uniqueAnalysis$pairwise[[iPairwise]]
      doTheseGenotypes <- if (doPsib != "all")
        rowSums(thisPairwise$match$flags == 0) == 0
      else
        rep(TRUE, nrow(thisPairwise$match$flags))

      if (sum(doTheseGenotypes) > 0) {
        for (iGenotype in which(doTheseGenotypes)) {
          thisGenotype <- thisPairwise$match$multilocus[iGenotype, ]
          psib <- 1
          for (iLocus in 1:max(multilocusMap)) {
            thisLocus <- thisGenotype[which(multilocusMap == iLocus)]
            if (all(thisPairwise$match$flags[iGenotype, which(multilocusMap == iLocus)] == 1)) {
              if (length(thisLocus) == 1) {
                alleleA <- thisLocus
                pA <- clusterAnalysisAlleleFreq$loci[[iLocus]]$alleleFreq[alleleA]
                psib <- psib * (1 + (2 * pA) + (pA * pA)) / 4
              } else {
                alleleA <- thisLocus[1]
                alleleB <- thisLocus[2]
                if (alleleA == alleleB) {
                  pA <- clusterAnalysisAlleleFreq$loci[[iLocus]]$alleleFreq[alleleA]
                  psib <- psib * (1 + (2 * pA) + (pA * pA)) / 4
                } else {
                  pA <- clusterAnalysisAlleleFreq$loci[[iLocus]]$alleleFreq[alleleA]
                  pB <- clusterAnalysisAlleleFreq$loci[[iLocus]]$alleleFreq[alleleB]
                  psib <- psib * (1 + pA + pB + (2 * pA * pB)) / 4
                }
              }
            }
          }
          uniqueAnalysis$pairwise[[iPairwise]]$match$psib[iGenotype] <- psib
          uniqueAnalysis$pairwise[[iPairwise]]$match$psibNotCalculable <- sum(!doTheseGenotypes)
        }
      }
    }

    indexUnclassified <- amDatasetFocal$index[!(amDatasetFocal$index %in% unique(do.call(
      c, lapply(clusterAnalysisPairwise$pairwise, function(x)
        c(x$focal$index, x$match$index))
    )))]
    uniqueAnalysis$numUnclassified <- length(indexUnclassified)

    if (uniqueAnalysis$numUnclassified > 0) {
      uniqueAnalysis$unclassified <- amDatasetFocal
      unclassifiedDatasetFocal <- which(amDatasetFocal$index %in% indexUnclassified)
      uniqueAnalysis$unclassified$index <- uniqueAnalysis$unclassified$index[unclassifiedDatasetFocal]
      if (!is.null(uniqueAnalysis$unclassified$metaData)) {
        uniqueAnalysis$unclassified$metaData <- uniqueAnalysis$unclassified$metaData[unclassifiedDatasetFocal]
      }
      tmpMultilocus <- matrix(
        uniqueAnalysis$unclassified$multilocus[unclassifiedDatasetFocal, ],
        length(indexUnclassified),
        ncol(uniqueAnalysis$unclassified$multilocus),
        byrow = FALSE
      )
      dimnames(tmpMultilocus)[[2]] <- dimnames(uniqueAnalysis$unclassified$multilocus)[[2]]
      uniqueAnalysis$unclassified$multilocus <- tmpMultilocus
    }

    indexMatches <- do.call(c,
                            lapply(clusterAnalysisPairwise$pairwise, function(x)
                              x$match$index))
    indexMultipleMatches <- unique(indexMatches[duplicated(indexMatches)])
    uniqueAnalysis$numMultipleMatches <- length(indexMultipleMatches)

    for (i in 1:length(uniqueAnalysis$pairwise)) {
      theseMultipleMatches <- which(uniqueAnalysis$pairwise[[i]]$match$index %in% indexMultipleMatches)
      if (length(theseMultipleMatches) != 0) {
        uniqueAnalysis$pairwise[[i]]$match$rowFlag[theseMultipleMatches] <- "MULTIPLE_MATCH"
        uniqueAnalysis$pairwise[[i]]$focal$rowFlag <- "CHECK_UNIQUE"
      } else {
        uniqueAnalysis$pairwise[[i]]$focal$rowFlag <- "UNIQUE"
      }
    }

    if (length(indexMultipleMatches) > 0) {
      uniqueAnalysis$multipleMatches <- amDatasetFocal
      multipleMatchesDatasetFocal <- which(amDatasetFocal$index %in% indexMultipleMatches)
      uniqueAnalysis$multipleMatches$index <- uniqueAnalysis$multipleMatches$index[multipleMatchesDatasetFocal]
      if (!is.null(uniqueAnalysis$multipleMatches$metaData)) {
        uniqueAnalysis$multipleMatches$metaData <- uniqueAnalysis$multipleMatches$metaData[multipleMatchesDatasetFocal]
      }
      tmpMultilocus <- matrix(
        uniqueAnalysis$multipleMatches$multilocus[multipleMatchesDatasetFocal, ],
        length(indexMultipleMatches),
        ncol(uniqueAnalysis$multipleMatches$multilocus),
        byrow = FALSE
      )
      dimnames(tmpMultilocus)[[2]] <- dimnames(uniqueAnalysis$multipleMatches$multilocus)[[2]]
      uniqueAnalysis$multipleMatches$multilocus <- tmpMultilocus
    }

    uniqueAnalysis$unique <- clusterAnalysis$unique
    uniqueAnalysis$unique$psib <- unlist(lapply(uniqueAnalysis$pairwise, function(x)
      x$match$psib[1]))
    uniqueAnalysis$unique$rowFlag <- unlist(lapply(uniqueAnalysis$pairwise, function(x)
      x$focal$rowFlag))
    uniqueAnalysis$cutHeight <- cutHeight
    uniqueAnalysis$consensusMethod <- clusterAnalysis$consensusMethod
    uniqueAnalysis$alleleMismatch <- alleleMismatch
    uniqueAnalysis$doPsib <- doPsib
    uniqueAnalysis$alleleFreq <- clusterAnalysisAlleleFreq

    return(uniqueAnalysis)
  }


#### amUniqueProfile() ####
amUniqueProfile <-
  function(amDatasetFocal,
           multilocusMap = NULL,
           alleleMismatch = NULL,
           matchThreshold = NULL,
           cutHeight = NULL,
           guessOptimum = TRUE,
           doPlot = TRUE,
           consensusMethod = 1,
           verbose = TRUE) {
    if (!inherits(amDatasetFocal, "amDataset"))
      stop("allelematch: amDatasetFocal must be an object of class \"amDataset\"",
           call. = FALSE)

    if (is.null(multilocusMap)) {
      if ((ncol(amDatasetFocal$multilocus) %% 2) != 0) {
        stop(
          "allelematch: there are an odd number of genotype columns in amDatasetFocal; Please specify multilocusMap manually",
          call. = FALSE
        )
      } else if (verbose) {
        cat(
          "allelematch: assuming genotype columns are in pairs, representing ",
          ncol(amDatasetFocal$multilocus) / 2,
          " loci\n",
          sep = ""
        )
      }
      multilocusMap <- rep(1:(ncol(amDatasetFocal$multilocus) / 2), each = 2)
    } else if (length(multilocusMap) != ncol(amDatasetFocal$multilocus)) {
      stop("allelematch: multilocusMap is the wrong length", call. = FALSE)
    } else if (sum(table(multilocusMap) > 2) > 0) {
      stop(
        "allelematch: multilocusMap indicates that a locus occurs in three or more columns",
        call. = FALSE
      )
    }
    multilocusMap <- as.integer(as.factor(multilocusMap))

    if (sum(!(c(
      is.null(alleleMismatch),
      is.null(matchThreshold),
      is.null(cutHeight)
    ))) > 1) {
      stop(
        "allelematch: please specify alleleMismatch OR matchThreshold OR cutHeight",
        call. = FALSE
      )
    }

    if (sum(!(c(
      is.null(alleleMismatch),
      is.null(matchThreshold),
      is.null(cutHeight)
    ))) == 0) {
      alleleMismatch <- seq(0, floor(length(multilocusMap)) * 0.4, 1)
      matchThreshold <- 1 - (alleleMismatch / length(multilocusMap))
      cutHeight <- 1 - matchThreshold
      profileType <- "alleleMismatch"
    } else {
      if (length(c(alleleMismatch, matchThreshold, cutHeight)) < 2)
        stop(
          "allelematch: please provide a range of parameter values for alleleMismatch OR matchThreshold OR cutHeight. e.g. alleleMismatch = c(0,1,2,3,4,5,6,7,8)",
          call. = FALSE
        )

      if (!is.null(alleleMismatch)) {
        if (length(alleleMismatch) <= 2)
          stop(
            "allelematch: alleleMismatch must be a vector containing a sequence of three or more values",
            call. = FALSE
          )
        matchThreshold <- 1 - (alleleMismatch / length(multilocusMap))
        cutHeight <- 1 - matchThreshold
        profileType <- "alleleMismatch"
      } else if (!is.null(cutHeight)) {
        if ((any(cutHeight < 0)) ||
            (any(cutHeight > 1)))
          stop("allelematch: cutHeight must be between 0 and 1", call. = FALSE)
        if (length(cutHeight) <= 2)
          stop("allelematch: cutHeight must contain three or more values",
               call. = FALSE)
        matchThreshold <- 1 - cutHeight
        alleleMismatch <- round((1 - matchThreshold) * length(multilocusMap), 2)
        profileType <- "cutHeight"
      } else {
        if ((any(matchThreshold < 0)) ||
            (any(matchThreshold > 1)))
          stop("allelematch: matchThreshold must be between 0 and 1",
               call. = FALSE)
        if (length(matchThreshold) <= 2)
          stop("allelematch: matchThreshold must contain three or more values",
               call. = FALSE)
        cutHeight <- 1 - matchThreshold
        alleleMismatch <- round((1 - matchThreshold) * length(multilocusMap), 2)
        profileType <- "matchThreshold"
      }
    }

    if (verbose)
      cat(
        "allelematch: running amUnique() at ",
        length(matchThreshold),
        " different values of ",
        profileType,
        "\n",
        sep = ""
      )

    profileResults <- data.frame(
      matchThreshold = NA,
      cutHeight = NA,
      alleleMismatch = NA,
      samples = NA,
      unique = NA,
      unclassified = NA,
      multipleMatch = NA,
      guessOptimum = NA
    )

    for (i in 1:length(matchThreshold)) {
      if (verbose) {
        cat(
          "allelematch: ",
          i,
          " of ",
          length(matchThreshold),
          " (matchThreshold=",
          round(matchThreshold[i], 2),
          ", cutHeight=",
          round(cutHeight[i], 2),
          ", alleleMismatch=",
          alleleMismatch[i],
          ")\n",
          sep = ""
        )
      }
      profileResults[i, "matchThreshold"] <- matchThreshold[i]
      profileResults[i, "cutHeight"] <- cutHeight[i]
      profileResults[i, "alleleMismatch"] <- alleleMismatch[i]
      amUniqueResult <- amUnique(
        amDatasetFocal,
        matchThreshold = matchThreshold[i],
        multilocusMap = multilocusMap,
        verbose = FALSE
      )
      profileResults[i, "unclassified"] <- amUniqueResult$numUnclassified
      profileResults[i, "multipleMatch"] <- amUniqueResult$numMultipleMatch
      profileResults[i, "unique"] <- amUniqueResult$focalDatasetN
      profileResults[i, "samples"] <- amUniqueResult$comparisonDatasetN
      profileResults[i, "guessOptimum"] <- NA
    }
    profileResults <- profileResults[order(matchThreshold, decreasing = TRUE), ]

    if (guessOptimum) {
      LeftTrimZero <- function(trimString) {
        oldString <- trimString
        repeat {
          newString <- sub("^0", "", oldString)
          if (newString == oldString)
            break
          oldString <- newString
        }
        return(newString)
      }

      i <- 0
      repeat {
        i <- i + 1
        minimum <- which.min(profileResults$multipleMatch[i:nrow(profileResults)])[1] + (i - 1)
        slopeAtMinimum <- if (minimum == 1)
          0
        else
          profileResults$multipleMatch[minimum] - profileResults$multipleMatch[minimum - 1]
        if ((slopeAtMinimum < 0) ||
            (minimum == nrow(profileResults)))
          break
      }
      secondMinimum <- minimum

      checkMultipleMatch <- LeftTrimZero(paste0(profileResults$multipleMatch, collapse = ""))
      leadingZeroes <- nchar(paste0(profileResults$multipleMatch, collapse = "")) - nchar(checkMultipleMatch)

      if (leadingZeroes == length(profileResults$multipleMatch)) {
        profileMorphology <- "ZeroFlat"
      } else if (sum(profileResults$multipleMatch[(leadingZeroes + 1):length(profileResults$multipleMatch)] == 0) == 0) {
        profileMorphology <- if (secondMinimum == length(profileResults$multipleMatch))
          "NoSecondMinimum"
        else
          "NonZeroSecondMinimum"
      } else {
        profileMorphology <- "ZeroSecondMinimum"
      }

      if (profileMorphology %in% c("ZeroFlat", "NoSecondMinimum")) {
        guessOptimum <- which.min(abs(sapply(2:(length(profileResults$unique) - 1), function(x)
          (profileResults$unique[x - 1] - profileResults$unique[x + 1]) / 2))) + 1
      } else {
        guessOptimum <- secondMinimum
      }

      profileResults$missingDataLoad <- round(
        sum(amDatasetFocal$multilocus == amDatasetFocal$missingCode) / length(amDatasetFocal$multilocus),
        3
      )
      profileResults$allelicDiversity <- round(mean(unlist(
        lapply(amUniqueResult$alleleFreq$loci, function(x)
          x$numAlleles)
      )), 1)
      profileResults$guessOptimum <- rep(FALSE, nrow(profileResults))
      profileResults$guessOptimum[guessOptimum] <- TRUE
      profileResults$guessMorphology <- profileMorphology

      if (verbose) {
        caution_text <- if (profileMorphology %in% c("NonZeroSecondMinimum", "NoSecondMinimum"))
          "allelematch: Use extra caution. Detection of optimal parameter is more error prone with this morphology\n"
        else
          ""
        cat(
          "allelematch: missing data load for input dataset is ",
          profileResults$missingDataLoad[1],
          "\n",
          "allelematch: allelic diversity for input dataset is ",
          profileResults$allelicDiversity[1],
          "\n",
          "allelematch: Best guess for optimal parameter at alleleMismatch=",
          profileResults$alleleMismatch[guessOptimum],
          " OR matchThreshold=",
          profileResults$matchThreshold[guessOptimum],
          " OR cutHeight=",
          profileResults$cutHeight[guessOptimum],
          "\n",
          "allelematch: Best guess for unique profile morphology: ",
          profileMorphology,
          "\n",
          caution_text,
          sep = ""
        )
      }
    }

    if (doPlot) {
      graphics::layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2), 1, 10))
      graphics::par(mar = c(5.1, 4.1, 2, 2))

      graphics::plot.default(
        c(min(profileResults[, profileType]), max(profileResults[, profileType])),
        c(0, max(profileResults[, c("unique", "unclassified", "multipleMatch")])),
        type = "n",
        axes = TRUE,
        xlab = profileType,
        ylab = "Count",
        cex = 2
      )
      graphics::points(profileResults[, profileType],
                       profileResults[, "unique"],
                       pch = 19,
                       col = "red")
      graphics::lines(
        profileResults[, profileType],
        profileResults[, "unique"],
        lwd = 2,
        lty = "solid",
        col = "red"
      )
      graphics::points(profileResults[, profileType],
                       profileResults[, "multipleMatch"],
                       pch = 22,
                       col = "black")
      graphics::lines(
        profileResults[, profileType],
        profileResults[, "multipleMatch"],
        lwd = 1,
        lty = "solid",
        col = "black"
      )
      graphics::points(profileResults[, profileType],
                       profileResults[, "unclassified"],
                       pch = 24,
                       col = "black")
      graphics::lines(
        profileResults[, profileType],
        profileResults[, "unclassified"],
        lwd = 1,
        lty = "dotted",
        col = "black"
      )

      if (guessOptimum) {
        graphics::arrows(
          profileResults[, profileType][which(profileResults$guessOptimum)],
          max(profileResults[, c("unique", "unclassified", "multipleMatch")]) * 0.3,
          profileResults[, profileType][which(profileResults$guessOptimum)],
          0,
          length = 0.15,
          angle = 20,
          lwd = 2
        )
      }

      graphics::par(mar = c(0, 0, 0, 1))
      graphics::plot.default(
        c(0, 100),
        c(0, 100),
        type = "n",
        axes = FALSE,
        ylab = "",
        xlab = ""
      )
      graphics::text(0,
                     100,
                     "allelematch",
                     cex = 2,
                     adj = c(0, 0.5))
      graphics::text(0,
                     97,
                     "amUniqueProfile()",
                     cex = 1,
                     adj = c(0, 0.5))
      graphics::text(
        0,
        88,
        paste0("missingDataLoad=", profileResults$missingDataLoad[1]),
        cex = 1,
        adj = c(0, 0.5)
      )
      graphics::text(
        0,
        86,
        paste0("allelicDiversity=", profileResults$allelicDiversity[1]),
        cex = 1,
        adj = c(0, 0.5)
      )
      graphics::legend(
        x = 0,
        y = 50,
        legend = c("unique", "multipleMatch", "unclassified"),
        lwd = c(2, 1, 1),
        lty = c("solid", "solid", "dotted"),
        col = c("red", "black", "black"),
        pch = c(19, 22, 24)
      )

      if (guessOptimum) {
        graphics::text(0,
                       35,
                       "Best guess for optimum:",
                       cex = 1,
                       adj = c(0, 0.5))
        graphics::text(
          0,
          33,
          paste0(profileType, "=", signif(profileResults[, profileType][which(profileResults$guessOptimum)], 2)),
          cex = 1,
          adj = c(0, 0.5)
        )
        graphics::text(0,
                       25,
                       "Profile morphology:",
                       cex = 1,
                       adj = c(0, 0.5))
        graphics::text(0,
                       23,
                       profileMorphology,
                       cex = 1,
                       adj = c(0, 0.5))
        if (profileMorphology %in% c("NonZeroSecondMinimum", "NoSecondMinimum")) {
          graphics::text(
            0,
            21,
            "Caution with optimum",
            col = "red",
            cex = 1,
            adj = c(0, 0.5)
          )
        }
      }
    }
    return(profileResults)
  }


#### summary.amUnique() ####
summary.amUnique <- function(object,
                             html = NULL,
                             csv = NULL,
                             ...) {
  if (!inherits(object, "amUnique"))
    stop("allelematch: this function requires an \"amUnique\" object",
         call. = FALSE)

  if (!is.null(html)) {
    if (is.logical(html) &&
        (html == TRUE))
      amHTML.amUnique(object)
    else if (is.logical(html) &&
             (html == FALSE))
      html <- NULL
    else
      amHTML.amUnique(object, htmlFile = html)
  }
  if (!is.null(csv))
    amCSV.amUnique(object, csvFile = csv, ...)

  if (is.null(html) && is.null(csv)) {
    cat(
      "allelematch: Console summary is not available for \"amUnique\" objects. Please use summary.amUnique(x, html=TRUE) or summary.amUnique(x, csv=\"file.csv\") options\n"
    )
  }
}


#### amCSSForHTML() ####
amCSSForHTML <- function() {
  return(
    "
        html { height: 100%; }
        body { background-color: inherit; color: inherit; font-family: Verdana; font-size: xx-small; margin: 0; height: 100%; }
        a:active, a:link, a:visited { color: #CC0000; }
        .amMismatchAllele { background-color: #CC0000; color: white; font-weight: bold; }
        .amMissingAllele { background-color: #FFCCCC; }
        .amInterpolatedAllele { background-color: blue; color: white; }
        .amGrid { border-collapse: separate; }
        .amGridContent { padding: 0; border: 1px solid #7EACB1; }
        .amGridUpperPanel, .amGridLowerPanel { padding: 3px; background-color: #F4FAFB; color: #2A769D; font-family: Verdana; font-size: xx-small; }
        .amGridUpperPanel { border-top: 0px; border-bottom: 1px solid; border-color: #7EACB1; }
        .amGridMiddlePanel { border: 0; }
        .amGridLowerPanel { border-top: 1px solid; border-bottom: 0px; border-color: #C2D4DA; }
        .amGridUpperPanel td, .amGridLowerPanel td { color: #2A769D; font-family: Verdana; font-size: xx-small; }
        .amTable { border: 0; border-spacing: 0; border-collapse: collapse; empty-cells: show; width: 100%; font-family: Verdana; font-size: xx-small; }
        .amTableSeparate { border-collapse: separate; }
        .amTable td { padding: 3px; border-bottom: 1px solid; border-right: 1px solid; border-color: #C2D4DA; white-space:nowrap; }
        .amTable .amTableHeader, .amTable .amTableHeader td { background-color: #B7D8DC; color: #000000; border-bottom: 1px solid; border-right: 1px solid; border-color: #7EACB1; vertical-align: top; white-space:nowrap; }
        .amPointer { cursor: pointer; }
        .amTableHeaderBtn { width: 100%; font-family: Verdana; font-size: xx-small; }
        .amTableHeader .amTableHeaderBtn td { background: transparent; padding: 0; border: 0; white-space: nowrap; }
        .amTableSelectRow { background-color: #FFFF66; color: #000000; }
    "
  )
}


#### amHTML.amPairwise() ####
amHTML.amPairwise <-
  function(x,
           htmlFile = NULL,
           htmlCSS = amCSSForHTML()) {
    if (!inherits(x, "amPairwise"))
      stop("allelematch: this function requires an \"amPairwise\" object",
           call. = FALSE)

    if (is.null(htmlFile)) {
      htmlFilePath <- tempdir()
      htmlFile <- file.path(htmlFilePath, "amPairwise.html")
      usingTmpFile <- TRUE
    } else {
      usingTmpFile <- FALSE
    }

    fileID <- file(htmlFile, "w")

    cat(
      "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n",
      "<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'><head>\n",
      "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></meta>\n",
      "<title>allelematch: ",
      class(x),
      "() output</title>\n",
      paste0(
        c("<style type=\"text/css\">", htmlCSS, "</style>"),
        collapse = "\n"
      ),
      "\n",
      "</head><body>\n",
      "<div style=\"margin-left:5%; margin-right:5%\">\n",
      "<br><table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">\n",
      "<div class=\"amGridUpperPanel\">\n",
      "<div style=\"font-size:x-small;\"><table>\n",
      "<tr><td style=\"width:400px;\"><span style=\"font-size:20px;\">allelematch<br>pairwise analysis</span><br><br></td></tr>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )

    headerHTML <- matrix("", 6, 2)
    headerHTML[1, ] <- c("\nfocal dataset N=", x$focalDatasetN)
    headerHTML[2, ] <- if (x$focalIsComparison)
      c("focal dataset compared against itself", "")
    else
      c("comparison dataset N=", x$comparisonDatasetN)
    headerHTML[3, ] <- c("missing data represented by: ", x$missingCode)
    headerHTML[4, ] <- c(
      "alleleMismatch (m-hat; maximum number of mismatching alleles): ",
      round(x$alleleMismatch, 2)
    )
    headerHTML[5, ] <- c(
      "matchThreshold (s-hat; lowest matching score returned): ",
      round(x$matchThreshold, 3)
    )
    headerHTML[6, ] <- c("summary generated: ", date())

    cat(
      paste0(
        "<tr><td><b>",
        headerHTML[, 1],
        "</b><em>",
        headerHTML[, 2],
        "</em></td></tr>\n",
        collapse = ""
      ),
      "</table></div></div></td></tr></table><br><br>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )

    y <- x$pairwise
    for (i in 1:length(y)) {
      headerNames <- if (!is.null(y[[i]]$focal$metaData))
        c("", "", dimnames(y[[i]]$focal$multilocus)[[2]], "Score")
      else
        c("", dimnames(y[[i]]$focal$multilocus)[[2]], "Score")
      headerTDs <- paste0(
        "<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
        headerNames,
        "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n",
        collapse = ""
      )

      cat(
        "<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">",
        "<div class=\"amGridUpperPanel\">",
        paste0(
          "<span style=\"font-size:medium;\">(",
          i,
          " of ",
          length(y),
          ")</span></div>"
        ),
        "<div class=\"amGridMiddlePanel\">\n<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n<tr class=\"amTableHeader\">\n",
        headerTDs,
        "</tr>\n<tr>\n",
        file = fileID,
        sep = "",
        append = TRUE
      )

      focalRow <- if (!is.null(y[[i]]$focal$metaData))
        c(y[[i]]$focal$index,
          y[[i]]$focal$metaData,
          y[[i]]$focal$multilocus,
          "FOCAL")
      else
        c(y[[i]]$focal$index, y[[i]]$focal$multilocus, "FOCAL")

      cat(
        paste0(
          "<td class=\"amTableSelectRow\"><div>",
          focalRow,
          "</div></td>\n",
          collapse = ""
        ),
        "</tr>\n",
        file = fileID,
        sep = "",
        append = TRUE
      )

      num_matches <- nrow(y[[i]]$match$multilocus)
      if (num_matches > 0) {
        match_rows_html <- character(num_matches)
        for (k in 1:num_matches) {
          if (!is.null(y[[i]]$focal$metaData)) {
            matchRow <- c(
              y[[i]]$match$index[k],
              y[[i]]$match$metaData[k],
              y[[i]]$match$multilocus[k, ],
              y[[i]]$match$score[k]
            )
            matchFlags <- c(1, 1, y[[i]]$match$flags[k, ], 1)
          } else {
            matchRow <- c(y[[i]]$match$index[k],
                          y[[i]]$match$multilocus[k, ],
                          y[[i]]$match$score[k])
            matchFlags <- c(1, y[[i]]$match$flags[k, ], 1)
          }

          cells <- character(length(matchRow))
          cells[matchFlags == 1] <- paste0("<td><div>", matchRow[matchFlags == 1], "</div></td>")
          cells[matchFlags == 0] <- paste0("<td><div class=\"amMismatchAllele\">",
                                           matchRow[matchFlags == 0],
                                           "</div></td>")
          cells[matchFlags == 2] <- paste0("<td><div class=\"amMissingAllele\">",
                                           matchRow[matchFlags == 2],
                                           "</div></td>")
          match_rows_html[k] <- paste0("<tr>", paste0(cells, collapse = ""), "</tr>\n")
        }
        cat(paste0(match_rows_html, collapse = ""),
            file = fileID,
            append = TRUE)
      }

      match_summary_text <- if ((y[[i]]$match$perfect == 0) &&
                                (y[[i]]$match$partial == 0))
        "<span>No matches found.</span>"
      else
        paste0(
          "<span>",
          y[[i]]$match$perfect,
          " perfect matches found. ",
          y[[i]]$match$partial,
          " partial matches found.</span>"
        )

      cat(
        "</table></div><div class=\"amGridLowerPanel\">",
        match_summary_text,
        "</div></td></tr></table><br><br>",
        file = fileID,
        sep = "",
        append = TRUE
      )
    }

    cat(
      "<span style=\"font-size:x-small;\">Generated by allelematch: an R package<br></span><span style=\"font-size:x-small;\">To reference this analysis please use citation(\"allelematch\")<br><br></span></div></body></html>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )
    close(fileID)

    if (usingTmpFile) {
      cat("Opening HTML file (", htmlFile, ") in default browser...\n", sep = "")

      # Only open the HTML if the session is interactive AND the amregtest suite
      # (or the user) hasn't explicitly disabled it via this environment variable.
      if (interactive() && Sys.getenv("ALLELEMATCH_SKIP_HTML") == "") {
        utils::browseURL(htmlFile)
      }

      oldTmpFiles <- Sys.glob(file.path(htmlFilePath, "am*.htm*"))
      if (length(oldTmpFiles) > 0)
        suppressWarnings(file.remove(oldTmpFiles))
    }
  }


#### amHTML.amCluster() ####
amHTML.amCluster <-
  function(x,
           htmlFile = NULL,
           htmlCSS = amCSSForHTML()) {
    if (!inherits(x, "amCluster"))
      stop("allelematch: this function requires an \"amCluster\" object",
           call. = FALSE)

    if (is.null(htmlFile)) {
      htmlFilePath <- tempdir()
      htmlFile <- file.path(htmlFilePath, "amCluster.html")
      usingTmpFile <- TRUE
    } else {
      usingTmpFile <- FALSE
    }

    fileID <- file(htmlFile, "w")

    cat(
      "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n",
      "<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'><head>\n",
      "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></meta>\n",
      "<title>allelematch: ",
      class(x),
      "() output</title>\n",
      paste0(
        c("<style type=\"text/css\">", htmlCSS, "</style>"),
        collapse = "\n"
      ),
      "\n",
      "</head><body>\n<div style=\"margin-left:5%; margin-right:5%\">\n",
      "<br><table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">\n<div class=\"amGridUpperPanel\">\n<div style=\"font-size:x-small;\"><table>\n",
      "<tr><td style=\"width:400px;\"><span style=\"font-size:20px;\">allelematch<br>cluster analysis</span><br><br></td></tr>\n",
      file = fileID,
      sep = ""
    )

    headerHTML <- matrix("", 12, 2)
    headerHTML[1, ] <- c("\nFocal dataset N=", x$focalDatasetN)
    headerHTML[2, ] <- c("unique N=", nrow(x$unique$multilocus))
    headerHTML[3, ] <- c("unique (consensus) N=", length(x$cluster))
    headerHTML[4, ] <- c("unique (singletons) N=", length(x$singletons))
    headerHTML[5, ] <- c("missing data represented by: ", x$missingCode)
    headerHTML[6, ] <- c("missing data matching method: ", x$missingMethod)
    headerHTML[7, ] <- c("clustered genotypes consensus method: ", x$consensusMethod)
    headerHTML[8, ] <- c("hierarchical clustering method: ", x$clusterMethod)
    headerHTML[9, ] <- c(
      "cutHeight (d-hat; dynamic tree cutting height): ",
      format(x$cutHeight, scientific = FALSE)
    )
    headerHTML[10, ] <- c("run until only singletons: ", x$runUntilSingletons)
    headerHTML[11, ] <- c("runs: ", x$totalRuns)
    headerHTML[12, ] <- c("summary generated: ", date())

    cat(
      paste0(
        "<tr><td><b>",
        headerHTML[, 1],
        "</b><em>",
        headerHTML[, 2],
        "</em></td></tr>\n",
        collapse = ""
      ),
      "</table></div></div></td></tr></table><br><br>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )

    y <- x$unique
    headerNames <- if (!is.null(y$metaData))
      c("", "", dimnames(y$multilocus)[[2]], "Type")
    else
      c("", dimnames(y$multilocus)[[2]], "Type")
    headerTDs <- paste0(
      "<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
      headerNames,
      "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n",
      collapse = ""
    )

    cat(
      "<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\"><div class=\"amGridUpperPanel\"><span style=\"font-size:medium;\">Unique genotypes</span></div>",
      "<div class=\"amGridMiddlePanel\">\n<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n<tr class=\"amTableHeader\">\n",
      headerTDs,
      "</tr>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )

    unique_rows_html <- character(nrow(y$multilocus))
    for (k in 1:nrow(y$multilocus)) {
      matchRow <- if (!is.null(y$metaData)) {
        c(
          paste0("<a href=\"#", y$index[k], "\">Jump</a>"),
          paste0("<a name=\"back_", y$index[k], "\">", y$index[k], "</a>"),
          y$metaData[k],
          y$multilocus[k, ],
          y$uniqueType[k]
        )
      } else {
        c(
          paste0("<a href=\"#", y$index[k], "\">Jump</a>"),
          paste0("<a name=\"back_", y$index[k], "\">", y$index[k], "</a>"),
          y$multilocus[k, ],
          y$uniqueType[k]
        )
      }
      unique_rows_html[k] <- paste0(
        "<tr onmouseout=\"this.style.cssText='background-color:none;';\" onmouseover=\"this.style.cssText='background-color:#E0FFFF';\">\n",
        paste0("<td><div>", matchRow, "</div></td>\n", collapse = ""),
        "</tr>\n"
      )
    }

    cat(
      paste0(unique_rows_html, collapse = ""),
      "</table></div><div class=\"amGridLowerPanel\"><span>There were ",
      nrow(y$multilocus),
      " unique genotypes identified using the parameters supplied.</span></div></td></tr></table><br><br>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )

    generate_group_html <- function(data_sub, type_label, cell_label) {
      if (length(data_sub) == 0)
        return("")
      out_blocks <- character(length(data_sub))

      for (i in 1:length(data_sub)) {
        headerNames <- if (!is.null(data_sub[[i]]$focal$metaData))
          c("", "", dimnames(data_sub[[i]]$focal$multilocus)[[2]], "Score")
        else
          c("", dimnames(data_sub[[i]]$focal$multilocus)[[2]], "Score")
        headerTDs <- paste0(
          "<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
          headerNames,
          "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n",
          collapse = ""
        )

        if (!is.null(data_sub[[i]]$focal$metaData)) {
          focalRow <- c(
            cell_label,
            data_sub[[i]]$focal$index,
            data_sub[[i]]$focal$metaData,
            data_sub[[i]]$focal$multilocus,
            ""
          )
          focalFlags <- c(1, 1, 1, data_sub[[i]]$focal$flags[1, ], 1)
        } else {
          focalRow <- c(cell_label,
                        data_sub[[i]]$focal$index,
                        data_sub[[i]]$focal$multilocus,
                        "")
          focalFlags <- c(1, 1, data_sub[[i]]$focal$flags[1, ], 1)
        }

        fcells <- character(length(focalRow))
        f_3 <- focalFlags == 3
        fcells[f_3] <- paste0(
          "<td class=\"amTableSelectRow\"><div class=\"amInterpolatedAllele\">",
          focalRow[f_3],
          "</div></td>"
        )
        fcells[!f_3] <- paste0("<td class=\"amTableSelectRow\"><div>",
                               focalRow[!f_3],
                               "</div></td>")
        f_row_text <- paste0("<tr>", paste0(fcells, collapse = ""), "</tr>\n")

        m_rows_text <- ""
        num_m <- nrow(data_sub[[i]]$match$multilocus)
        if (num_m > 0) {
          m_lines <- character(num_m)
          for (k in 1:num_m) {
            if (!is.null(data_sub[[i]]$focal$metaData)) {
              matchRow <- c(
                if (cell_label == "CONSENSUS")
                  "MEMBER"
                else
                  "CLOSEST",
                data_sub[[i]]$match$index[k],
                data_sub[[i]]$match$metaData[k],
                data_sub[[i]]$match$multilocus[k, ],
                data_sub[[i]]$match$score[k]
              )
              matchFlags <- c(1, 1, 1, data_sub[[i]]$match$flags[k, ], 1)
            } else {
              matchRow <- c(
                if (cell_label == "CONSENSUS")
                  "MEMBER"
                else
                  "CLOSEST",
                data_sub[[i]]$match$index[k],
                data_sub[[i]]$match$multilocus[k, ],
                data_sub[[i]]$match$score[k]
              )
              matchFlags <- c(1, 1, data_sub[[i]]$match$flags[k, ], 1)
            }

            cells <- character(length(matchRow))
            cells[matchFlags == 1] <- paste0("<td><div>", matchRow[matchFlags == 1], "</div></td>")
            cells[matchFlags == 0] <- paste0("<td><div class=\"amMismatchAllele\">",
                                             matchRow[matchFlags == 0],
                                             "</div></td>")
            cells[matchFlags == 2] <- paste0("<td><div class=\"amMissingAllele\">",
                                             matchRow[matchFlags == 2],
                                             "</div></td>")
            m_lines[k] <- paste0("<tr>", paste0(cells, collapse = ""), "</tr>\n")
          }
          m_rows_text <- paste0(m_lines, collapse = "")
        }

        footer_info <- if (cell_label == "CONSENSUS") {
          if ((data_sub[[i]]$match$perfect == 0) &&
              (data_sub[[i]]$match$partial == 0))
            "<span>No matches found.</span>"
          else
            paste0(
              "<span>",
              data_sub[[i]]$match$perfect,
              " perfect matches found. ",
              data_sub[[i]]$match$partial,
              " partial matches found.</span>"
            )
        } else {
          "<span>Closest match to singleton shown for diagnostic purposes</span>"
        }

        out_blocks[i] <- paste0(
          "<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\"><div class=\"amGridUpperPanel\"><a name=\"",
          data_sub[[i]]$focal$index,
          "\"><span style=\"font-size:medium;\">(",
          type_label,
          " ",
          i,
          " of ",
          length(data_sub),
          ")</span></a><br><a href=\"#back_",
          data_sub[[i]]$focal$index,
          "\">Jump up</a><br></div><div class=\"amGridMiddlePanel\">\n<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n<tr class=\"amTableHeader\">\n",
          headerTDs,
          "</tr>\n<tr>\n",
          f_row_text,
          m_rows_text,
          "</table></div><div class=\"amGridLowerPanel\">",
          footer_info,
          "</div></td></tr></table><br><br>\n"
        )
      }
      return(paste0(out_blocks, collapse = ""))
    }

    cat(
      generate_group_html(x$cluster, "Cluster", "CONSENSUS"),
      file = fileID,
      append = TRUE
    )
    cat(
      generate_group_html(x$singletons, "Singleton", "SINGLETON"),
      file = fileID,
      append = TRUE
    )

    cat(
      "<span style=\"font-size:x-small;\">Generated by allelematch: an R package<br></span><span style=\"font-size:x-small;\">To reference this analysis please use citation(\"allelematch\")<br><br></span></div>\n</body></html>",
      file = fileID,
      sep = "",
      append = TRUE
    )
    close(fileID)

    if (usingTmpFile) {
      cat("Opening HTML file (", htmlFile, ") in default browser...\n", sep = "")

      # Only open the HTML if the session is interactive AND the amregtest suite
      # (or the user) hasn't explicitly disabled it via this environment variable.
      if (interactive() && Sys.getenv("ALLELEMATCH_SKIP_HTML") == "") {
        utils::browseURL(htmlFile)
      }

      oldTmpFiles <- Sys.glob(file.path(htmlFilePath, "am*.htm*"))
      if (length(oldTmpFiles) > 0)
        suppressWarnings(file.remove(oldTmpFiles))
    }
  }


#### amHTML.amUnique() ####
amHTML.amUnique <-
  function(x,
           htmlFile = NULL,
           htmlCSS = amCSSForHTML()) {
    if (!inherits(x, "amUnique"))
      stop("allelematch: this function requires an \"amUnique\" object",
           call. = FALSE)

    if (is.null(htmlFile)) {
      htmlFilePath <- tempdir()
      htmlFile <- file.path(htmlFilePath, "amUnique.html")
      usingTmpFile <- TRUE
    } else {
      usingTmpFile <- FALSE
    }

    fileID <- file(htmlFile, "w")

    cat(
      "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Transitional//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd\">\n",
      "<html xmlns='http://www.w3.org/1999/xhtml' xml:lang='en' lang='en'><head>\n",
      "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></meta>\n",
      "<title>allelematch: ",
      class(x),
      "() output</title>\n",
      paste0(
        c("<style type=\"text/css\">", htmlCSS, "</style>"),
        collapse = "\n"
      ),
      "\n",
      "</head><body>\n<div style=\"margin-left:5%; margin-right:5%\">\n",
      "<br><table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\">\n<div class=\"amGridUpperPanel\">\n<div style=\"font-size:x-small;\"><table>\n",
      "<tr><td style=\"width:400px;\"><span style=\"font-size:20px;\">allelematch<br>unique analysis</span><br><br></td></tr>\n",
      file = fileID,
      sep = ""
    )

    if (inherits(x, "amUnique")) {
      headerHTML <- matrix("", 14, 2)
      headerHTML[1, ] <- c("\nunique N=", x$focalDatasetN)
      headerHTML[2, ] <- c("samples N=", x$comparisonDatasetN)
      headerHTML[3, ] <- c("<br>loci N=", length(x$alleleFreq$loci))
      headerHTML[4, ] <- c("locus names: ", paste(sapply(x$alleleFreq$loci, function(z)
        z$name), collapse = ", "))
      headerHTML[5, ] <- c("<br>missing data represented by: ", x$missingCode)
      headerHTML[6, ] <- c("clustered genotypes consensus method: ", x$consensusMethod)
      headerHTML[7, ] <- c(
        "Psib calculated for: ",
        ifelse(
          x$doPsib != "all",
          "samples with no mismatches",
          "all samples (mismatches treated as missing data)"
        )
      )
      headerHTML[8, ] <- c(
        "<br>alleleMismatch (m-hat; maximum number of mismatching alleles): ",
        round(x$alleleMismatch, 2)
      )
      headerHTML[9, ] <- c("cutHeight (d-hat; dynamic tree cutting height): ",
                           round(x$cutHeight, 3))
      headerHTML[10, ] <- c(
        "matchThreshold (s-hat; lowest matching score returned): ",
        round(x$matchThreshold, 3)
      )
      headerHTML[11, ] <- if (x$numUnclassified > 0)
        c(
          "<br>unclassified N=",
          paste0(
            "<span style=\"color:red;\">",
            x$numUnclassified,
            "</span>"
          )
        )
      else
        c("<br>unclassified N=", x$numUnclassified)
      headerHTML[12, ] <- c("multipleMatch N=", x$numMultipleMatches)
      headerHTML[13, ] <- c(
        "<br>Note: Unique genotypes are determined based on clustering of scores.",
        " Psib appears for reference purposes."
      )
      headerHTML[14, ] <- c("<br>summary generated: ", date())
    }

    cat(
      paste0(
        "<tr><td><b>",
        headerHTML[, 1],
        "</b><em>",
        headerHTML[, 2],
        "</em></td></tr>\n",
        collapse = ""
      ),
      "</table></div></div></td></tr></table><br><br>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )

    y <- x$unique
    headerNames <- if (!is.null(y$metaData))
      c("", "", "", dimnames(y$multilocus)[[2]], "Psib", "Type")
    else
      c("", "", dimnames(y$multilocus)[[2]], "Psib", "Type")
    headerTDs <- paste0(
      "<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
      headerNames,
      "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n",
      collapse = ""
    )

    cat(
      "<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\"><div class=\"amGridUpperPanel\"><span style=\"font-size:medium;\">Unique genotypes</span></div>",
      "<div class=\"amGridMiddlePanel\">\n<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n<tr class=\"amTableHeader\">\n",
      headerTDs,
      "</tr>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )

    unique_rows_html <- character(nrow(y$multilocus))
    for (k in 1:nrow(y$multilocus)) {
      psib <- signif(y$psib[k], 2)
      psib <- if (is.na(psib))
        "&nbsp; &nbsp; ---"
      else if (psib < 0.00001)
        "<0.00001"
      else
        format(psib, scientific = FALSE)

      matchRow <- if (!is.null(y$metaData)) {
        c(
          paste0("<a href=\"#", y$index[k], "\">Jump</a>"),
          y$index[k],
          y$metaData[k],
          y$multilocus[k, ],
          psib,
          y$rowFlag[k]
        )
      } else {
        c(
          paste0("<a href=\"#", y$index[k], "\">Jump</a>"),
          y$index[k],
          y$multilocus[k, ],
          psib,
          y$rowFlag[k]
        )
      }
      unique_rows_html[k] <- paste0(
        "<tr onmouseout=\"this.style.cssText='background-color:none;';\" onmouseover=\"this.style.cssText='background-color:#E0FFFF';\">\n",
        paste0("<td><div>", matchRow, "</div></td>\n", collapse = ""),
        "</tr>\n"
      )
    }

    cat(
      paste0(unique_rows_html, collapse = ""),
      "</table></div><div class=\"amGridLowerPanel\"><span>There were ",
      nrow(y$multilocus),
      " unique genotypes identified using the parameters supplied.</span></div></td></tr></table><br><br>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )

    if (!is.null(x$unclassified)) {
      y <- x$unclassified
      headerNames <- if (!is.null(y$metaData))
        c("", "", "", dimnames(y$multilocus)[[2]], "Type")
      else
        c("", "", dimnames(y$multilocus)[[2]], "Type")
      headerTDs <- paste0(
        "<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
        headerNames,
        "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n",
        collapse = ""
      )

      cat(
        "<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\"><div class=\"amGridUpperPanel\"><span style=\"font-size:medium; color:red; font-weight:bold;\">Unclassified samples</span></div><div class=\"amGridMiddlePanel\">\n<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n<tr class=\"amTableHeader\">\n",
        headerTDs,
        "</tr>\n",
        file = fileID,
        sep = "",
        append = TRUE
      )

      unclass_rows_html <- character(nrow(y$multilocus))
      for (k in 1:nrow(y$multilocus)) {
        matchRow <- if (!is.null(y$metaData)) {
          c(
            "<span style=\"color:red; font-weight:bold;\">!!!</span>",
            y$index[k],
            y$metaData[k],
            y$multilocus[k, ],
            "UNCLASSIFIED"
          )
        } else {
          c(
            "<span style=\"color:red; font-weight:bold;\">!!!</span>",
            y$index[k],
            y$multilocus[k, ],
            "UNCLASSIFIED"
          )
        }
        unclass_rows_html[k] <- paste0(
          "<tr onmouseout=\"this.style.cssText='background-color:none;';\" onmouseover=\"this.style.cssText='background-color:#E0FFFF';\">\n",
          paste0("<td><div>", matchRow, "</div></td>\n", collapse = ""),
          "</tr>\n"
        )
      }

      cat(
        paste0(unclass_rows_html, collapse = ""),
        "</table></div><div class=\"amGridLowerPanel\"><span>There were ",
        nrow(y$multilocus),
        " samples that were not classified.</span></div></td></tr></table><br><br>\n",
        file = fileID,
        sep = "",
        append = TRUE
      )
    }

    y <- x$pairwise
    for (i in 1:length(y)) {
      headerNames <- if (!is.null(y[[i]]$focal$metaData))
        c("",
          "",
          "",
          dimnames(y[[i]]$focal$multilocus)[[2]],
          "Psib",
          "Score",
          "Type")
      else
        c("",
          "",
          dimnames(y[[i]]$focal$multilocus)[[2]],
          "Psib",
          "Score",
          "Type")
      headerTDs <- paste0(
        "<td class=\"amPointer\"><table cellspacing=\"0\" class=\"amTableHeaderBtn\"><tr><td>",
        headerNames,
        "</td><td style=\"width:10px;\">&nbsp;</td></tr></table></td>\n",
        collapse = ""
      )

      cat(
        "<table cellspacing=\"0\" class=\"amGrid\"><tr><td class=\"amGridContent\"><div class=\"amGridUpperPanel\"><a name=\"",
        y[[i]]$focal$index,
        "\"><span style=\"font-size:medium;\">Unique genotype (",
        i,
        " of ",
        length(y),
        ")</span></a></div><div class=\"amGridMiddlePanel\">\n<table cellspacing=\"0\" class=\"amTable amTableSeparate\">\n<tr class=\"amTableHeader\">\n",
        headerTDs,
        "</tr>\n<tr>\n",
        file = fileID,
        sep = "",
        append = TRUE
      )

      focalRow <- if (!is.null(y[[i]]$focal$metaData))
        c(
          "",
          y[[i]]$focal$index,
          y[[i]]$focal$metaData,
          y[[i]]$focal$multilocus,
          "",
          "",
          y[[i]]$focal$rowFlag
        )
      else
        c("",
          y[[i]]$focal$index,
          y[[i]]$focal$multilocus,
          "",
          "",
          y[[i]]$focal$rowFlag)

      cat(
        paste0(
          "<td class=\"amTableSelectRow\"><div>",
          focalRow,
          "</div></td>\n",
          collapse = ""
        ),
        "</tr>\n",
        file = fileID,
        sep = "",
        append = TRUE
      )

      num_matches <- nrow(y[[i]]$match$multilocus)
      if (num_matches > 0) {
        match_rows_html <- character(num_matches)
        for (k in 1:num_matches) {
          psib <- signif(y[[i]]$match$psib[k], 2)
          psib <- if (is.na(psib))
            "&nbsp; &nbsp; ---"
          else if (psib < 0.00001)
            "<0.00001"
          else
            format(psib, scientific = FALSE)

          is_multimatch <- y[[i]]$match$rowFlag[k] == "MULTIPLE_MATCH"
          lead_flag <- if (is_multimatch)
            "<span style=\"color:red; font-weight:bold;\">!!!</span>"
          else
            ""
          type_lbl  <- if (is_multimatch)
            "MULTIPLE_MATCH"
          else
            "MATCH"

          if (!is.null(y[[i]]$focal$metaData)) {
            matchRow <- c(
              lead_flag,
              y[[i]]$match$index[k],
              y[[i]]$match$metaData[k],
              y[[i]]$match$multilocus[k, ],
              psib,
              y[[i]]$match$score[k],
              type_lbl
            )
            matchFlags <- c(1, 1, 1, y[[i]]$match$flags[k, ], 1, 1, 1)
          } else {
            matchRow <- c(
              lead_flag,
              y[[i]]$match$index[k],
              y[[i]]$match$multilocus[k, ],
              psib,
              y[[i]]$match$score[k],
              type_lbl
            )
            matchFlags <- c(1, 1, y[[i]]$match$flags[k, ], 1, 1, 1)
          }

          cells <- character(length(matchRow))
          cells[matchFlags == 1] <- paste0("<td><div>", matchRow[matchFlags == 1], "</div></td>")
          cells[matchFlags == 0] <- paste0("<td><div class=\"amMismatchAllele\">",
                                           matchRow[matchFlags == 0],
                                           "</div></td>")
          cells[matchFlags == 2] <- paste0("<td><div class=\"amMissingAllele\">",
                                           matchRow[matchFlags == 2],
                                           "</div></td>")
          match_rows_html[k] <- paste0("<tr>", paste0(cells, collapse = ""), "</tr>\n")
        }
        cat(paste0(match_rows_html, collapse = ""),
            file = fileID,
            append = TRUE)
      }

      cat(
        "</table></div><div class=\"amGridLowerPanel\">",
        file = fileID,
        append = TRUE
      )

      if ((y[[i]]$match$perfect == 0) &&
          (y[[i]]$match$partial == 0)) {
        cat("<span>No matches found.</span>",
            file = fileID,
            append = TRUE)
      } else {
        multipleMatchMsg <- if (any(y[[i]]$match$rowFlag == "MULTIPLE_MATCH"))
          "<br><span style=\"color:red\">multipleMatch samples present</span>"
        else
          ""
        doPsibMsg <- if (x$doPsib != "all")
          "<br>Psib calculated for samples with no mismatches."
        else
          "<br>Psib calculated for all samples."
        cat(
          paste0(
            "<span>Unique genotype compared against ",
            x$comparisonDatasetN,
            " samples, returning score>=",
            round(x$matchThreshold, 3),
            multipleMatchMsg,
            doPsibMsg,
            "</span>"
          ),
          file = fileID,
          sep = "",
          append = TRUE
        )
      }
      cat("</div>\n</td></tr></table><br><br>\n",
          file = fileID,
          append = TRUE)
    }

    cat(
      "<span style=\"font-size:x-small;\">Generated by allelematch: an R package<br></span><span style=\"font-size:x-small;\">To reference this analysis please use citation(\"allelematch\")<br><br></span></div>\n</body></html>\n",
      file = fileID,
      sep = "",
      append = TRUE
    )
    close(fileID)

    if (usingTmpFile) {
      cat("Opening HTML file (", htmlFile, ") in default browser...\n", sep = "")

      # Only open the HTML if the session is interactive AND the amregtest suite
      # (or the user) hasn't explicitly disabled it via this environment variable.
      if (interactive() && Sys.getenv("ALLELEMATCH_SKIP_HTML") == "") {
        utils::browseURL(htmlFile)
      }

      oldTmpFiles <- Sys.glob(file.path(htmlFilePath, "am*.htm*"))
      if (length(oldTmpFiles) > 0)
        suppressWarnings(file.remove(oldTmpFiles))
    }
  }


#### amCSV.amPairwise() ####
amCSV.amPairwise <- function(x, csvFile) {
  if (!inherits(x, "amPairwise"))
    stop("allelematch: this function requires an \"amPairwise\" object",
         call. = FALSE)

  y <- x$pairwise
  has_metadata <- !is.null(y[[1]]$match$metaData)
  match_list <- vector("list", length(y))
  match_threshold_str <- as.character(x$matchThreshold)

  for (i in 1:length(y)) {
    match_group_count <- nrow(y[[i]]$match$multilocus)
    if (has_metadata) {
      match_list[[i]] <- cbind(
        matchGroup = rep(i, match_group_count),
        nMatchGroup = rep(match_group_count, match_group_count),
        focalIndex = rep(y[[i]]$focal$index, match_group_count),
        comparisonIndex = y[[i]]$match$index,
        matchThreshold = rep(match_threshold_str, match_group_count),
        score = y[[i]]$match$score,
        comparisonMetaData = y[[i]]$match$metaData,
        y[[i]]$match$multilocus
      )
    } else {
      match_list[[i]] <- cbind(
        matchGroup = rep(i, match_group_count),
        nMatchGroup = rep(match_group_count, match_group_count),
        focalIndex = rep(y[[i]]$focal$index, match_group_count),
        comparisonIndex = y[[i]]$match$index,
        matchThreshold = rep(match_threshold_str, match_group_count),
        score = y[[i]]$match$score,
        y[[i]]$match$multilocus
      )
    }
  }

  csvTable <- do.call(rbind, match_list)
  if (is.null(dim(csvTable)))
    csvTable <- matrix(csvTable,
                       nrow = 1,
                       dimnames = list(NULL, names(csvTable)))
  utils::write.csv(csvTable, file = csvFile, row.names = FALSE)
}


#### amCSV.amCluster() ####
amCSV.amCluster <- function(x, csvFile) {
  if (!inherits(x, "amCluster"))
    stop("allelematch: this function requires an \"amCluster\" object",
         call. = FALSE)
  y <- x$unique
  csvTable <- if (!is.null(y$metaData))
    data.frame(
      index = as.character(y$index),
      metaData = as.character(y$metaData),
      y$multilocus,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  else
    data.frame(
      index = as.character(y$index),
      y$multilocus,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  utils::write.csv(csvTable, file = csvFile, row.names = FALSE)
}


#### amCSV.amUnique() ####
amCSV.amUnique <- function(x, csvFile, uniqueOnly = FALSE) {
  if (!inherits(x, "amUnique"))
    stop("allelematch: this function requires an \"amUnique\" object",
         call. = FALSE)

  if (uniqueOnly) {
    class(x) <- "amCluster"
    amCSV.amCluster(x, csvFile)
  } else {
    y <- x$pairwise
    has_metadata <- !is.null(y[[1]]$match$metaData)
    out_list <- vector("list", length(y) * 2 + 1)
    list_idx <- 1

    if (x$numUnclassified > 0) {
      u <- x$unclassified
      u_cnt <- nrow(u$multilocus)
      u_df <- data.frame(
        uniqueGroup = rep("", u_cnt),
        rowType = rep("UNCLASSIFIED", u_cnt),
        uniqueIndex = as.character(u$index),
        matchIndex = rep("", u_cnt),
        nUniqueGroup = rep("", u_cnt),
        alleleMismatch = rep(as.character(x$alleleMismatch), u_cnt),
        matchThreshold = rep(as.character(x$matchThreshold), u_cnt),
        cutHeight = rep(as.character(x$cutHeight), u_cnt),
        Psib = rep("", u_cnt),
        score = rep("", u_cnt),
        stringsAsFactors = FALSE
      )
      if (!is.null(u$metaData))
        u_df$metaData <- as.character(u$metaData)
      out_list[[list_idx]] <- cbind(u_df, u$multilocus)
      list_idx <- list_idx + 1
    }

    mismatch_str <- as.character(x$alleleMismatch)
    threshold_str <- as.character(x$matchThreshold)
    cutheight_str <- as.character(x$cutHeight)

    for (i in 1:length(y)) {
      rowFlag <- y[[i]]$match$rowFlag
      rowFlag[rowFlag == ""] <- "MATCH"
      Psib <- y[[i]]$match$psib
      Psib[is.na(Psib)] <- ""
      match_cnt <- nrow(y[[i]]$match$multilocus)

      focal_df <- data.frame(
        uniqueGroup = i,
        rowType = as.character(y[[i]]$focal$rowFlag),
        uniqueIndex = as.character(y[[i]]$focal$index),
        matchIndex = as.character(y[[i]]$focal$index),
        nUniqueGroup = match_cnt,
        alleleMismatch = mismatch_str,
        matchThreshold = threshold_str,
        cutHeight = cutheight_str,
        Psib = as.character(y[[i]]$match$psib[1]),
        score = "",
        stringsAsFactors = FALSE
      )
      if (has_metadata)
        focal_df$metaData <- as.character(y[[i]]$focal$metaData)
      out_list[[list_idx]] <- cbind(focal_df, y[[i]]$focal$multilocus)
      list_idx <- list_idx + 1

      match_mask <- y[[i]]$match$index != y[[i]]$focal$index
      sub_cnt <- sum(match_mask)
      if (sub_cnt > 0) {
        match_df <- data.frame(
          uniqueGroup = rep(i, sub_cnt),
          rowType = as.character(rowFlag[match_mask]),
          uniqueIndex = rep(as.character(y[[i]]$focal$index), sub_cnt),
          matchIndex = as.character(y[[i]]$match$index[match_mask]),
          nUniqueGroup = rep(match_cnt, sub_cnt),
          alleleMismatch = rep(mismatch_str, sub_cnt),
          matchThreshold = rep(threshold_str, sub_cnt),
          cutHeight = rep(cutheight_str, sub_cnt),
          Psib = as.character(Psib[match_mask]),
          score = as.character(y[[i]]$match$score[match_mask]),
          stringsAsFactors = FALSE
        )
        if (has_metadata)
          match_df$metaData <- as.character(y[[i]]$match$metaData[match_mask])
        out_list[[list_idx]] <- cbind(match_df, y[[i]]$match$multilocus[match_mask, , drop = FALSE])
        list_idx <- list_idx + 1
      }
    }

    out_list <- out_list[1:(list_idx - 1)]
    csvTable <- do.call(rbind, out_list)
    utils::write.csv(csvTable, file = csvFile, row.names = FALSE)
  }
}


#### amPreCheck() ####
amPreCheck <- function(amDatasetFocal) {
  if (!inherits(amDatasetFocal, "amDataset")) {
    stop(
      "allelematch: amPreCheck: amDatasetFocal must be an object of class 'amDataset'; use amDataset() first",
      call. = FALSE
    )
  }

  genotypes <- amDatasetFocal$multilocus
  missing_code <- amDatasetFocal$missingCode

  missing_per_ind <- rowSums(genotypes == missing_code) / ncol(genotypes)
  max_ind_missing <- max(missing_per_ind) * 100
  mean_ind_missing <- mean(missing_per_ind) * 100

  missing_per_locus <- colSums(genotypes == missing_code) / nrow(genotypes)
  max_locus_missing <- max(missing_per_locus) * 100

  # 1. Identify pairs (assuming columns are strictly ordered a, b, a, b...)
  n_cols <- ncol(genotypes)

  # Determine pairs: integer division plus remainder for odd columns
  # e.g., 21 columns -> 10 pairs + 1 single = 11 loci total
  n_loci <- (n_cols %/% 2) + (n_cols %% 2)

  # 2. Create a locus-level validity matrix (True if at least one allele is called)
  # This collapses the 2-column diploid pair into 1 logical locus
  locus_validity <- matrix(NA, nrow = nrow(genotypes), ncol = n_loci)

  for (i in 1:(n_cols %/% 2)) {
    col_a <- (genotypes[, 2 * i - 1] != missing_code)
    col_b <- (genotypes[, 2 * i]     != missing_code)
    locus_validity[, i] <- (col_a | col_b)
  }

  # If there is an odd dangling column, handle it as a standalone locus
  if (n_cols %% 2 != 0) {
    locus_validity[, n_loci] <- (genotypes[, n_cols] != missing_code)
  }

  # Calculate cross-product matrix based on the corrected locus count
  shared_loci_matrix <- tcrossprod(locus_validity)
  diag(shared_loci_matrix) <- NA
  min_overlap <- min(shared_loci_matrix, na.rm = TRUE)

  if (min_overlap == 0) {
    verdict_str <- paste0(
      "  [CRITICAL WARNING] At least two individuals exhibit ZERO overlapping loci.\n",
      "  Running amUniqueProfile() or amCluster() may result in an atomic sorting failure or the creation of statistically unstable clusters.\n",
      "  Action Required: Verify the clustering output carefully or pre-filter loci or individuals with high missing data rates before running allelematch.\n"
    )
  } else if (max_locus_missing > 50 || mean_ind_missing > 35) {
    verdict_str <- paste0(
      "  [CAUTION] Missing data load is heavy. High parameter settings may produce unstable clusters.\n",
      "  Action Suggested: Consider conservative quality control screening.\n"
    )
  } else {
    verdict_str <- "  [SAFE] Dataset looks complete enough for clustering.\n"
  }

  cat(
    "\n========================================================\n",
    "          allelematch Data Pre-Screening Quality Report   \n",
    "========================================================\n",
    "Locus:\n",
    "  - Maximum missing data found in a single locus: ",
    round(max_locus_missing, 1),
    "%\n\n",
    "Individual Sample:\n",
    "  - Average missing data load across all individuals: ",
    round(mean_ind_missing, 1),
    "%\n",
    "  - Worst individual missing data load profile:       ",
    round(max_ind_missing, 1),
    "%\n\n",
    "Overlap Check:\n",
    "  - Minimum shared loci between any pair of samples:  ",
    min_overlap,
    " out of ",
    n_loci,
    "\n\n",
    "Verdict:\n",
    verdict_str,
    "========================================================\n\n",
    sep = ""
  )

  invisible(
    list(
      maxLocusMissing = max_locus_missing,
      meanIndMissing = mean_ind_missing,
      minOverlap = min_overlap
    )
  )
}
