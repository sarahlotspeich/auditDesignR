% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_bcc.R
\name{sample_bcc}
\alias{sample_bcc}
\title{Sample unvalidated balanced case-control (BCC*) audit based on Phase I variables.}
\usage{
sample_bcc(dat, phI, phII, sample_on)
}
\arguments{
\item{dat}{Dataframe or matrix containing columns \code{sample_on}.}

\item{phI}{Phase I sample size.}

\item{phII}{Phase II sample size.}

\item{sample_on}{Columns with the Phase I variables (should be categorical) used for sampling strata (can be name or numeric index). Currently, sampling on up to 3 binary variables can be accommodated.}
}
\value{
A vector of length \code{phI} with validation indicators V = 1 if selected for Phase II and V = 0 otherwise.
}
\description{
Sample unvalidated balanced case-control (BCC*) audit based on Phase I variables.
}
