% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_fn.R
\name{setchecknames}
\alias{setchecknames}
\title{Set check.names attribute of omicsData object}
\usage{
setchecknames(omicsData, set_to = TRUE)
}
\arguments{
\item{omicsData}{an object of the class 'pepData', 'proData', 'metabData', 'lipidData', or 'nmrData', usually created by \code{\link{as.pepData}}, \code{\link{as.proData}}, \code{\link{as.metabData}}, \code{\link{as.lipidData}}, or \code{\link{as.nmrData}}, respectively.}

\item{set_to}{logical indicating what to set check.names attribute to. Defaults to TRUE.}
}
\value{
updated omicsData object with check.names attribute
}
\description{
This function sets the check.names attribute of an omicsData object
}
\examples{
dontrun{
setchecknames(omicsData, set_to = TRUE)
}
}
