% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/05_basis.R
\docType{class}
\name{fourier+-class}
\alias{fourier+-class}
\title{Fourier Basis Functions}
\description{
The "\code{fourier+}" class represents Fourier basis functions.
}
\section{Slots}{

\describe{
\item{\code{range}}{A length-two numeric vector to define the interval on which
basis functions can be evaluated. Default is \code{c(0, 1)}.}

\item{\code{nbasis}}{The total number of basis functions, including the ones that need
to be dropped. Default is 3.}

\item{\code{dropind}}{Indices of basis functions that need to be dropped. Default is
\code{numeric(0)} (no basis will be dropped).}

\item{\code{ncoef}}{The actual number of basis functions to represent functional data,
satisfying the condition \code{ncoef == nbasis - length(dropind)}.}

\item{\code{period}}{The length of the cycle that Fourier basis functions repeat
             themselves. Default is 1.}
}}
\section{Extends}{

Class \code{\link[=basis+-class]{basis+}}
}

