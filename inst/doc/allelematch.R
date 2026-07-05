## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----example1-load------------------------------------------------------------
# library(allelematch)
# 
# data(amExample1)

## ----example1-amDataset-amPreCheck--------------------------------------------
# example1 <- amDataset(
#   amExample1,
#   indexColumn = "sampleId",
#   ignoreColumn = "knownIndividual",
#   missingCode = "-99"
# )
# 
# amPreCheck(example1)

## ----example1-amUniqueProfile-------------------------------------------------
# amUniqueProfile(example1, doPlot = TRUE)

## ----example1-amUnique, comment = ""------------------------------------------
# uniqueExample1 <- amUnique(example1, alleleMismatch = 2)

## ----example1-summary-amUnique-htmlSummary------------------------------------
# ## Save to disk
# summary(uniqueExample1, html = "example1_1.html")
# 
# ## View in default browser
# summary(uniqueExample1, html = TRUE)

## ----example1-summary-amUnique-csvSummary-------------------------------------
# summary(uniqueExample1, csv = "example1_1.csv")

## ----example1-summary-amUnique-csvSummaryUniqueOnly---------------------------
# summary(uniqueExample1, csv = "example1_2.csv", uniqueOnly = TRUE)

## ----example1-amDataset-amUnique-htmlSummary----------------------------------
# example1chk <- amDataset(
#   amExample1,
#   indexColumn = "sampleId",
#   metaDataColumn = "knownIndividual",
#   missingCode = "-99"
# )
# 
# uniqueExample1chk <- amUnique(example1chk, alleleMismatch = 2)
# 
# summary(uniqueExample1chk, html = "example1_2.html")

## ----example2-amDataset-amPreCheck--------------------------------------------
# data(amExample2)
# 
# example2 <- amDataset(
#   amExample2,
#   indexColumn = "sampleId",
#   metaDataColumn = "knownIndividual",
#   missingCode = "-99"
# )
# 
# amPreCheck(example2)

## ----example2-amUniqueProfile-------------------------------------------------
# amUniqueProfile(example2, doPlot = TRUE)

## ----example2-amUnique-htmlSummary--------------------------------------------
# uniqueExample2 <- amUnique(example2, alleleMismatch = 3)
# 
# summary(uniqueExample2, html = "example2_1.html")

## ----example2-amUnique-psib-htmlSummary---------------------------------------
# uniqueExample2 <- amUnique(example2, alleleMismatch = 3, doPsib = "all")
# 
# summary(uniqueExample2, html = "example2_2.html")

## ----example3-amDataset-amPreCheck--------------------------------------------
# data(amExample3)
# 
# example3 <- amDataset(
#   amExample3,
#   indexColumn = "sampleId",
#   metaDataColumn = "knownIndividual",
#   missingCode = "-99"
# )
# 
# amPreCheck(example3)

## ----example3-amUniqueProfile-------------------------------------------------
# amUniqueProfile(example3, doPlot = TRUE)

## ----example3-amUnique-htmlSummary--------------------------------------------
# uniqueExample3 <- amUnique(example3, alleleMismatch = 6)
# 
# summary(uniqueExample3, html = "example3_1.html")

## ----example3-amPairwise------------------------------------------------------
# unclassifiedExample3 <- amPairwise(
#   uniqueExample3$unclassified,
#   uniqueExample3$unique,
#   alleleMismatch = 7
# )

## ----example3-htmlSummary-----------------------------------------------------
# summary(unclassifiedExample3, html = "example3_2.html")

## ----example3-amPairwise-htmlSummary------------------------------------------
# multipleMatchExample3 <- amPairwise(
#   uniqueExample3$multipleMatch,
#   uniqueExample3$unique,
#   alleleMismatch = 6
# )
# 
# summary(multipleMatchExample3, html = "example3_3.html")

## ----example3-csvSummary------------------------------------------------------
# summary(uniqueExample3, csv = "example3_1.csv")

## ----example4-amDataset-amPreCheck--------------------------------------------
# data(amExample4)
# 
# example4 <- amDataset(
#   amExample4,
#   indexColumn = "sampleId",
#   metaDataColumn = "knownIndividual",
#   missingCode = "-99"
# )
# 
# amPreCheck(example4)

## ----example4-amUniqueProfile-------------------------------------------------
# amUniqueProfile(example4, doPlot = TRUE)

## ----example4-1-amUnique_htmlSummary------------------------------------------
# uniqueExample4 <- amUnique(example4, alleleMismatch = 1)
# 
# summary(uniqueExample4, html = "example4_1.html")

## ----example4-2-amUnique-htmlSummary------------------------------------------
# uniqueExample4ballpark <- amUnique(example4, alleleMismatch = 6)
# 
# summary(uniqueExample4ballpark, html = "example4_2.html")

## ----example4-3-amUnique-htmlSummary------------------------------------------
# uniqueExample4high <- amUnique(example4, alleleMismatch = 8)
# 
# summary(uniqueExample4high, html = "example4_3.html")

## ----example5-amDataset-amPreCheck--------------------------------------------
# data(amExample5)
# 
# # Inspect the unique levels of metadata and identifiers
# head(levels(amExample5$samplingData))
# head(levels(amExample5$sampleId))
# head(levels(amExample5$sex))
# 
# example5 <- amDataset(
#   amExample5,
#   indexColumn = "sampleId",
#   metaDataColumn = "samplingData",
#   missingCode = "-99"
# )
# 
# amPreCheck(example5)

## ----example5-names-----------------------------------------------------------
# names(amExample5)

## ----example5-map1------------------------------------------------------------
# example5map <- c(
#   "sex", "LOC1", "LOC1", "LOC2", "LOC2", "LOC3", "LOC3",
#   "LOC4", "LOC4", "LOC5", "LOC5", "LOC6", "LOC6", "LOC7",
#   "LOC7", "LOC8", "LOC8", "LOC9", "LOC9", "LOC10", "LOC10"
# )

## ----example5-map2------------------------------------------------------------
# example5map <- c(1, rep(2:11, each = 2))

## ----example5-amUniqueProfile-------------------------------------------------
# amUniqueProfile(example5, multilocusMap = example5map, doPlot = TRUE)

## ----example5-1-amUnique-htmlSummary------------------------------------------
# uniqueExample5 <- amUnique(
#   example5,
#   multilocusMap = example5map,
#   alleleMismatch = 3
# )
# 
# summary(uniqueExample5, html = "example5_1.html")

## ----example5-2-amPairwise-htmlSummary----------------------------------------
# unclassifiedExample5 <- amPairwise(
#   uniqueExample5$unclassified,
#   uniqueExample5$unique,
#   alleleMismatch = 4
# )
# 
# summary(unclassifiedExample5, html = "example5_2.html")

## ----example5-3-amPairwise-htmlSummary----------------------------------------
# multipleMatchExample5 <- amPairwise(
#   uniqueExample5$multipleMatch,
#   uniqueExample5$unique,
#   alleleMismatch = 3
# )
# 
# summary(multipleMatchExample5, html = "example5_3.html")

## ----example5-csvSummary------------------------------------------------------
# summary(uniqueExample5, csv = "example5_1.csv")

