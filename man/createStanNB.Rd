% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{createStanNB}
\alias{createStanNB}
\title{Create a Negative Binomial .stan model}
\usage{
createStanNB(P = 2)
}
\arguments{
\item{P}{The highest power of WBS signal in regression}
}
\value{
The .stan model
}
\description{
This is an internal function used to create a Stan model.
}
\examples{
cat(createStanNB(2))
}
