% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TMB.R
\name{SR}
\alias{SR}
\title{Sequential reduction configuration}
\usage{
SR(x, discrete = FALSE)
}
\arguments{
\item{x}{Breaks defining the domain of integration}

\item{discrete}{Boolean defining integration wrt Lebesgue measure (\code{discrete=FALSE}) or counting measure \code{discrete=TRUE}.}
}
\description{
Helper function to specify an integration grid used by the
sequential reduction algorithm available through the argument
\code{integrate} to \code{MakeADFun}.
}
