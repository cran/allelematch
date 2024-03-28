\name{amUnique}
\alias{amUnique}
\alias{summary.amUnique}
\alias{amHTML.amUnique}
\alias{amCSV.amUnique}

\title{
	Identification of unique genotypes
}

\description{
	Identifies unique genotypes and generates analysis output in formatted text, HTML, or
	CSV. Samples are clustered and matched based on their dissimilarity score (see
	\code{\link{amMatrix}}). Also calculated is the match probability, Psib, which is the
	probability that a sample is a sibling of a unique genotype (and therefore not a
	replicate sample) given the allele frequencies in a population consisting of only the
	unique genotypes (Wilberg & Dreher, 2004).
}

\usage{
	amUnique(
		amDatasetFocal,
		multilocusMap = NULL,
		alleleMismatch = NULL,
		matchThreshold = NULL,
		cutHeight = NULL,
		doPsib = "missing",
		consensusMethod = 1,
		verbose = TRUE
		)

	amHTML.amUnique(
		x,
		htmlFile = NULL,
		htmlCSS = amCSSForHTML()
		)

	amCSV.amUnique(
		x,
		csvFile,
		uniqueOnly = FALSE
		)


  \method{summary}{amUnique}(
    object,
    html = NULL,
    csv = NULL,
    ...
    )
}

\arguments{
	\item{amDatasetFocal}{
		An \code{amDataset} object containing genotypes in which an unknown number of
		individuals are sampled multiple times
	}
	\item{multilocusMap}{
		Optional.\cr
		A vector of integers or strings giving the mappings onto loci for all genotype
		columns in amDatasetFocal.\cr
		When omitted, columns are assumed to be paired (i.e., diploid loci with alleles in
		adjacent columns).\cr
		See details.
	}
	\item{alleleMismatch}{
		Optional.\cr
		Maximum number of mismatching alleles which will be tolerated when identifying
		individuals; also known as m-hat parameter.\cr
		If specified, then \code{matchThreshold} and \code{cutHeight} should be
		omitted.\cr
		All three parameters are related.\cr
		See details.
	}
	\item{matchThreshold}{
		Optional.\cr
		Minimum dissimilarity score which constitutes a match when identifying
		individuals; also known as s-hat parameter.\cr
		If specified, then \code{alleleMismatch} and \code{cutHeight} should be omitted;
		all three parameters are related.\cr
		See details.
	}
	\item{cutHeight}{
		Optional.\cr
		The \code{cutHeight} parameter used in dynamic tree cutting by \code{amCluster};
		also known as d-hat parameter.\cr
		If specified, then \code{alleleMismatch} and \code{matchThreshold} should be
		omitted.\cr
		All three parameters are related.\cr
		See details.
	}
	\item{doPsib}{
		String specifying how match probability should be calculated.\cr
		See details.
	}
	\item{consensusMethod}{
		The method (an integer) used to determine the consensus multilocus genotype from a
		cluster of multilocus genotypes.\cr
		See \code{\link{amCluster}} for details.\cr
		Typically the default is adequate.
	}
	\item{verbose}{
		If \code{verbose = TRUE}, report the progress of the analysis to the console.\cr
		Useful with datasets consisting of thousands of samples where progress may be
		slow.
	}
	\item{object, x}{
		An \code{amUnique} object.
	}
	\item{htmlFile}{
		HTML filepath to create.\cr
		If \code{htmlFile = NULL}, a file is created in the operating system temporary
		directory and is then opened in the default browser.
	}
	\item{htmlCSS}{
		A string containing a valid cascading style sheet.\cr
		A default style sheet is provided in \code{\link{amCSSForHTML}}.\cr
		See \code{\link{amCSSForHTML}} for details of how to tweak this CSS.
	}
	\item{html}{
		If \code{html = TRUE}, the \code{summary.amUnique} method produces and loads an HTML file in the
		default browser.\cr
		\code{html} can also contain a path to a file where HTML output will be
		written.\cr
		Note that \code{summary.amUnique} does not produce formatted output for the
		console.
	}
	\item{csvFile, csv}{
		CSV filepath to create containing a representation of the \code{amUnique}
		analysis.
	}
	\item{uniqueOnly}{
		If \code{uniqueOnly = TRUE}, only the unique genotypes will be saved to a CSV,
		with no additional information associated with the analysis.
	}
	\item{...}{
		Additional arguments to \code{summary.amUnique}
	}
}

\details{
	Only one of \code{alleleMismatch}, \code{cutHeight}, \code{matchThreshold} can be
	specified, as the three parameters are related.\cr

	\code{alleleMismatch} is the most intuitive way to understand how the identification
	of unique genotypes proceeds. For example, a setting of \code{alleleMismatch = 4}
	implies that up to four alleles may be different for multiple samples to be
	representatives of the same individual. In practice, however, this value is only an
	approximation of the amount of mismatch that may be tolerated. This is because the
	clustering process used to identify unique genotypes, and the subsequent matching
	which identifies samples that match these unique genotypes is based on a dissimilarity
	metric or score (see \code{\link{amMatrix}}) that incorporates both allele mismatches
	and missing data. \code{alleleMismatch} is not used in analyses and is converted to
	this dissimilarity metric in the following manner: \code{cutHeight} which is parameter
	for \code{\link{amCluster}} and called from this function is \code{cutHeight =
	alleleMismatch/(number of allele columns)} and \code{matchThreshold} which is a
	parameter for \code{\link{amPairwise}} and also called from this function is
	\code{matchThreshold = 1 - cutHeight}.\cr

	Selecting the appropriate value for \code{alleleMismatch}, \code{cutHeight}, or
	\code{matchThreshold} is an important task. Use \code{\link{amUniqueProfile}} to
	assist in this process. Seethe Data S1 Supplementary documentation and tutorials (PDF)
	located at <doi:10.1111/j.1755-0998.2012.03137.x>\cr

	\code{doPsib = "missing"} is the default and specifies that match probability Psib
	should be calculated for samples that match unique genotypes and have no allele
	mismatches, but may differ by having missing data. \code{doPsib = "all"} specifies
	that Psib should be calculated for all samples that match unique genotypes. In this
	case, if allele mismatches occur, alleles are assumed to be missing at the mismatching
	loci.\cr

	\code{multilocusMap} is often not required, as amDataset objects will typically
	consist of paired columns of genotypes, where each pair is a separate locus. In cases
	where this is not the case (e.g., gender is in only one column), a map vector must be
	specified.\cr

	Example: \code{amDataset} consists of gender followed by 4 diploid loci in paired
	columns\cr
	\code{multilocusMap = c(1, 2, 2, 3, 3, 4, 4, 5, 5)}\cr
	or equally\cr
	\code{multilocusMap=c("GENDER", "LOC1", "LOC1", "LOC2", "LOC2", "LOC3", "LOC4",
	"LOC4")}

	For more information on selecting \code{consensusMethod} see \code{\link{amCluster}}.
	The default \code{consensusMethod = 1} is typically adequate.
}

\value{
	\code{amUnique} object or side effects: analysis summary written to an HTML file or to
	the console, or written to a CSV file.
}

\note{
	There is an additional side effect of \code{html = TRUE} (or of \code{htmlFile =
	NULL}). If required, there is a clean up of the operating system temporary directory
	where AlleleMatch temporary HTML files are stored. Files that match the pattern
	am*.html and are older 24 hours are deleted from this temporary directory.
}

\references{
	For a complete vignette, please access via the Data S1 Supplementary documentation and
	tutorials (PDF) located at <doi:10.1111/j.1755-0998.2012.03137.x>.\cr

	Wilberg MJ, Dreher BP (2004) GENECAP: a program for analysis of multilocus genotype
	data for non-invasive sampling and capture-recapture population estimation.
	\emph{Molecular Ecology Notes}, 4, 783-785.
}

\author{
	Paul Galpern (\email{pgalpern@gmail.com})
}

\seealso{
	\code{\link{amDataset}}, \code{\link{amMatrix}}, \code{\link{amPairwise}},
	\code{\link{amCluster}}, \code{\link{amUniqueProfile}}
}

\examples{
	\dontrun{
	data("amExample2")

	## Produce amDataset object
	myDataset <-
		amDataset(
			amExample2,
			missingCode = "-99",
			indexColumn = 1,
			ignoreColumn = 2
			)

	## Usage
	## Optimal alleleMismatch parameter previously found using amUniqueProfile()
	myUnique <-
		amUnique(
		myDataset,
		alleleMismatch = 3
		)

	## Display analysis as HTML in default browser
	summary.amUnique(
		myUnique,
		html = TRUE
		)

	## Save analysis to HTML file
	summary.amUnique(
		myUnique,
		html = "myUnique.htm"
		)

	## Save analysis to a CSV file
	summary.amUnique(
		myUnique,
		csv = "myUnique.csv"
		)

	## Save unique genotypes only to a CSV file
	summary.amUnique(
		myUnique,
		csv = "myUnique.csv",
		uniqueOnly = TRUE
		)

	## Data set with gender information
	data("amExample5")

	## Produce amDataset object
	myDataset2 <-
		amDataset(
			amExample5,
			missingCode = "-99",
			indexColumn = 1,
			metaDataColumn = 2
			)

	## Usage
	## Optimal alleleMismatch parameter previously found using amUniqueProfile()
	myUniqueProfile <-
		amUnique(
			myDataset2,
			multilocusMap = c(1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10,
			11, 11),
			alleleMismatch = 3
			)

	}

}
