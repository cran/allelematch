\name{amPairwise}
\alias{amPairwise}
\alias{summary.amPairwise}
\alias{amHTML.amPairwise}
\alias{amCSV.amPairwise}

\title{
	Pairwise matching of multilocus genotypes
}

\description{
	Functions to perform a pairwise matching analysis of a multilocus genotype dataset, and review
	the output in formatted text or HTML. For each genotype in the focal dataset all genotypes in
	the comparison genotype are returned that match at or above a threshold matching score. The
	matching score is also known as the s-hat criterion (see the supplementary documentation). This
	is determined using \code{\link{amMatrix}}.
}

\usage{
	amPairwise(
		amDatasetFocal,
		amDatasetComparison = amDatasetFocal,
		alleleMismatch = NULL,
		matchThreshold = NULL,
		missingMethod = 2
		)

	amHTML.amPairwise(
		x,
		htmlFile = NULL,
		htmlCSS = amCSSForHTML()
		)

	amCSV.amPairwise(
		x,
		csvFile
		)

  \method{summary}{amPairwise}(
    object,
    html = NULL,
    csv = NULL,
    ...
    )
}

\arguments{
	\item{amDatasetFocal}{
		An \code{amDataset} object containing focal genotypes.
	}
	\item{amDatasetComparison}{
		Optional.\cr
		An \code{amDataset} object containing comparison genotypes.\cr
		If not supplied, the focal dataset is also the comparison dataset (i.e., all focal dataset
		members are compared against one another).
	}
	\item{alleleMismatch}{
		Maximum number of mismatching alleles which will be tolerated when identifying individuals;
		also known as m-hat parameter.\cr
		If specified, then \code{matchThreshold} should be omitted.
	}
	\item{matchThreshold}{
		Return comparison genotypes that match with the focal genotype at or above this score or
		similarity; also known as s-hat parameter.
	}
	\item{missingMethod}{
		Method used to determine the similarity of multilocus genotypes when data is missing.\cr
		The default, (\code{missingMethod = 2}), is preferable in all cases.\cr
		See \code{\link{amMatrix}}.
	}
	\item{object, x}{
		An \code{amPairwise} object.
	}
	\item{htmlFile}{
		HTML filepath to create.\cr
		If \code{htmlFile = NULL}, a file is created in the operating system temporary directory and
		is then opened in the default browser.
	}
	\item{htmlCSS}{
		A string containing a valid cascading style sheet.\cr
		A default style sheet is provided in \code{\link{amCSSForHTML}}.\cr
		See \code{\link{amCSSForHTML}} for details of how to tweak this CSS.
	}
	\item{html}{
		If \code{html = NULL} or \code{html=FALSE}, formatted textual output is displayed on the
		console.\cr
		If \code{html = TRUE}, the \code{summary.amPairwise} method produces and loads an HTML file in the
		default browser.\cr
		\code{html} can also contain a path to a file where HTML output will be written.
	}
	\item{csvFile, csv}{
		CSV filepath to create containing giving a data frame representation of the pairwise
		matching results.
	}
	\item{...}{
		Additional arguments to \code{summary.amPairwise}.
	}
}

\details{
	Pairwise matching of genotypes is a useful means to assess data quality and inspect for
	genotyping errors.\cr

	\code{matchThreshold} represents the similarity between two multilocus genotypes and can be
	thought of as a percentage similarity (or a Hamming's distance between two vectors) that has
	been corrected where missing data is present, such that missing data represents neither a match
	nor a mismatch but a "partial" match. See \code{\link{amMatrix}} for more discussion of this
	metric.
}

\note{
	As \code{matchThreshold} is lowered, the size of the output increases rapidly. Typically
	analyses will not be very useful or manageable with thresholds below 0.7.\cr

	There is an additional side effect of \code{html = TRUE} (or of \code{htmlFile = NULL}). If
	required, there is a clean up of the operating system temporary directory where AlleleMatch
	temporary HTML files are stored. Files that match the pattern am*.html and are older 24 hours
	are deleted from this temporary directory.
}

\value{
	\code{amPairwise} object or side effects: analysis summary written to an HTML file or to the
	console, or written to a CSV file.
}

\references{
	For a complete vignette, please access via the Data S1 Supplementary documentation and tutorials
	(PDF) located at <doi:10.1111/j.1755-0998.2012.03137.x>.
}

\author{
	Paul Galpern (\email{pgalpern@gmail.com})
}

\seealso{
	\code{\link{amDataset}}, \code{\link{amMatrix}}, \code{\link{amUnique}}
}

\examples{
	\dontrun{
	data("amExample5")

	## Produce amDataset object
	myDataset <-
		amDataset(
			amExample5,
			missingCode = "-99",
			indexColumn = 1,
			metaDataColumn = 2,
			ignoreColumn = "gender"
			)

	## Typical usage
	myPairwise <-
		amPairwise(
			myDataset,
			alleleMismatch = 2
			)

	## Display analysis as HTML in default browser
	summary.amPairwise(
		myPairwise,
		html = TRUE
		)

	## Save analysis to HTML file
	summary.amPairwise(
		myPairwise,
		html = "myPairwise.htm"
		)

	## Save analysis to CSV file
	summary.amPairwise(
		myPairwise,
		csv = "myPairwise.csv"
		)

	## Display analysis as formatted text on the console
	summary.amPairwise(myPairwise)

	## Compare one dataset against a second
	## Both must have same number of allele columns
	## Here we create two datasets artificially from one for illustration purposes
	myDatasetA <-
		amDataset(
			amExample5[sample(nrow(amExample5))[1:25], ],
			missingCode = "-99",
			indexColumn = 1,
			ignoreColumn = 2
			)
	myDatasetB <-
		amDataset(
			amExample5[sample(nrow(amExample5))[1:100], ],
			missingCode = "-99",
			indexColumn = 1,
			ignoreColumn = 2
			)
	myPairwise2 <-
		amPairwise(
			myDatasetA,
			myDatasetB,
			alleleMismatch = 3
			)
	summary.amPairwise(
		myPairwise2,
		html = TRUE
		)
	}
}
