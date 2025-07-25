% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treedater0.R
\name{sampleYearsFromLabels}
\alias{sampleYearsFromLabels}
\title{Compute a vector of numeric sample times from labels in a sequence aligment or phylogeny}
\usage{
sampleYearsFromLabels(
  tips,
  dateFormat = "\%Y-\%m-\%d",
  delimiter = NULL,
  index = NULL,
  regex = NULL
)
}
\arguments{
\item{tips}{A character vector supplying the name of each sample}

\item{dateFormat}{The format of the sample date. See ?Date for more information}

\item{delimiter}{Character(s) which separate data in each label}

\item{index}{Integer position of the date string in each label with respect to *delimiter*}

\item{regex}{A regular expression for finding the date substring. Should not be used with *delimiter* or *index*}
}
\value{
Numeric vector with sample time in decimal format.
}
\description{
Compute a vector of numeric sample times from labels in a sequence aligment or phylogeny
}
\examples{
## A couple of labels for Ebola virus sequences:
sampleYearsFromLabels( c('EBOV|AA000000|EM104|SierraLeone_EM|2014-06-02'
                       , 'EBOV|AA000000|G3713|SierraLeone_G|2014-06-09')
, delimiter='|' )
## Equivalently:
sampleYearsFromLabels( c('EBOV|AA000000|EM104|SierraLeone_EM|2014-06-02'
                       , 'EBOV|AA000000|G3713|SierraLeone_G|2014-06-09')
 , regex='[0-9]+-[0-9]+-[0-9]+')

}
