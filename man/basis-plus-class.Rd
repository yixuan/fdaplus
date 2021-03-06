% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/05_basis.R
\docType{class}
\name{basis+-class}
\alias{basis+-class}
\title{Basis Functions for Functional Data}
\description{
The "\code{basis+}" class is the base class for all basis function objects,
including \code{\link[=bspline+-class]{bspline+}} and
\code{\link[=fourier+-class]{fourier+}}.
}
\section{Slots}{

\describe{
\item{\code{range}}{A length-two numeric vector to define the interval on which
basis functions can be evaluated.}

\item{\code{nbasis}}{The total number of basis functions, including the ones that need
to be dropped.}

\item{\code{dropind}}{Indices of basis functions that need to be dropped.}

\item{\code{ncoef}}{The actual number of basis functions to represent functional data,
            satisfying the condition \code{ncoef == nbasis - length(dropind)}.}
}}
\seealso{
\code{\link[=bspline+-class]{bspline+}},
\code{\link[=fourier+-class]{fourier+}}.
}

