% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/10_bspline.R
\name{basis_bs}
\alias{basis_bs}
\title{Creating B-spline Basis Functions}
\usage{
basis_bs(range = c(0, 1), nbasis = NULL, order = 4, breaks = NULL,
  dropind = NULL)
}
\arguments{
\item{range}{A length-two numeric vector to define the interval on which
basis functions can be evaluated. Default is \code{c(0, 1)}.}

\item{nbasis}{The total number of basis functions, including the ones that need
to be dropped. See section \strong{Details}.}

\item{order}{The order of the B-spline functions, which is the degree of
splines plus one. Default is 4, standing for cubic splines.}

\item{breaks}{A vector of break points (inner knots) to define the B-spline basis.
This does not include the end points, which are defined by
the parameter \code{range}. See section \strong{Details}.}

\item{dropind}{Indices of basis functions that need to be dropped. Default is
               \code{NULL}, meaning no basis will be dropped.}
}
\value{
A \code{\link[=bspline+-class]{bspline+}} object representing the
        basis functions.
}
\description{
This function constructs a \code{\link[=bspline+-class]{bspline+}} object
that represents a series of B-spline basis functions.
}
\details{
\code{nbasis}, \code{order} and \code{breaks} are related to each other,
satisfying the condition \code{nbasis == order + length(breaks)}. Hence when
any one of \code{nbasis} and \code{breaks} is properly specified,
the other one can be computed accordingly. For details,

\itemize{
  \item If both \code{nbasis} and \code{breaks} are \code{NULL} (the default),
        \code{nbasis} will be set to \code{order} and there is no inner knot.
  \item If \code{nbasis} is \code{NULL} and \code{breaks} is given, then
        \code{nbasis} will be calculated using the equation above.
  \item If \code{nbasis} is given and \code{breaks} unspecified, a series
        of equally spaced inner knots will be calculated.
  \item If both \code{nbasis} and \code{breaks} are provided, then the program
        will check their validity.
}
}
\author{
Yixuan Qiu <\url{http://statr.me/}>
}

