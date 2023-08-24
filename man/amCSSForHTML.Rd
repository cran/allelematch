\name{amCSSForHTML}
\alias{amCSSForHTML}

\title{
	Produce cascading style sheet (CSS) for HTML
}

\description{
	Returns a string containing the CSS code for embedding in HTML output by
	\code{\link{amHTML.amPairwise}} and \code{\link{amHTML.amUnique}}.
}

\usage{
	amCSSForHTML()
}

\details{
	This function is used internally. It can also be used as a basis to tweak the CSS
code if different output formatting and colour-coding are desired. See examples.
}

\value{
	A string containing a complete cascading style sheet.
}

\references{
	For more information and for a complete vignette, please access via the Data S1	Supplementary
	documentation and tutorials (PDF) located at <doi:10.1111/j.1755-0998.2012.03137.x>.
}

\author{
	Paul Galpern (\email{pgalpern@gmail.com})
}

\seealso{
	\code{\link{amHTML.amPairwise}}, \code{\link{amHTML.amUnique}}
}

\examples{
	\dontrun{
	data("amExample5")
	
	## Produce CSS file for editing
	cat(amCSSForHTML(), file = "myCSS.css")
	
	## Edit myCSS.css manually using text editor, then reimport as follows
	myCSS <- 
		paste(
			readLines("myCSS.css"), 
			collapse = "\n"
			)
	
	## Perform an allelematch unique analysis
	myUnique <- 
		amUnique(
			amDataset(
				amExample5, 
				missingCode = "-99", 
				indexColumn = 1, 
				metaDataColumn = 2, 
				ignoreColumn = "gender"
				),
			alleleMismatch = 3
			)

	## Produce HTML output using tweaked CSS
	amHTML.amUnique(
		myUnique, 
		html = "myUnique.htm", 
		htmlCSS = myCSS
		)
	}
}