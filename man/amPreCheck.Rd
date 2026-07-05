\name{amPreCheck}
\alias{amPreCheck}
\title{Pre-Screen an amDataset Object for Severe Missing Data Loads}
\description{
  Analyzes missing data footprints and pairwise locus overlaps within an amDataset object to anticipate mathematical failures during subsequent clustering.
}
\usage{
  amPreCheck(amDatasetFocal)
}
\arguments{
  \item{amDatasetFocal}{An object of class \code{amDataset}.}
}
\value{
  Prints a quality control metric dashboard to the console. Also returns an invisible list containing summary percentages.
}
\author{Todd Cross}
\examples{
  \dontrun{
    data(amExample1)
    techreps <- amDataset(amExample1, missingCode = "-99")
    amPreCheck(techreps)
  }
}
