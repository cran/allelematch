\name{amExampleData}
\alias{amExampleData}
\alias{amExample1}
\alias{amExample2}
\alias{amExample3}
\alias{amExample4}
\alias{amExample5}
\docType{data}
\title{
Data sets to support the tutorials in the supplementary documentation and examples in the manual.
}
\description{
\code{amExample1} High quality data set\cr
\code{amExample2} Good quality data set\cr
\code{amExample3} Marginal quality data set\cr
\code{amExample4} Low quality data set\cr
\code{amExample5} Wildlife data set\cr

Data sets 1 to 4 are simulated.  Because the data are simulated,
the individual from which the samples are derived is known.  This is given
in the column knownIndividual, and permits an assessment of how effective
the software has been.

Data set 5 is a real dataset from a wildlife population that has been non-invasively
sampled.  The meta-data has been fabricated.  It represents the output from allele calling software.


}
\usage{data(amExample1)
data(amExample2)
data(amExample3)
data(amExample4)
data(amExample5)}
\format{
  Data frames with differing numbers of samples in rows, and alleles in
  columns.  Missing data is represented as "-99".
  }
\details{
Note how in amExample5 a single marker "Gender" has been combined with diploid loci.\cr\cr
Also note how in all data sets in diploid loci the lower length allele always comes first.
This pattern is typical in allele calling software.
}