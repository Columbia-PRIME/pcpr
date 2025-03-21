% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparsity.R
\name{sparsity}
\alias{sparsity}
\title{Estimate sparsity of given matrix}
\usage{
sparsity(D, thresh = 1e-04)
}
\arguments{
\item{D}{The input data matrix.}

\item{thresh}{(Optional) A numeric threshold \verb{>= 0} used to determine if an
entry in \code{D} is "practically zero". If the absolute value of an entry is
below \code{thresh}, then it is judged to be "practically zero". By default,
\code{thresh = 1e-04}.}
}
\value{
The sparsity of \code{D}, measured as the percentage of entries in \code{D}
that are "practically zero".
}
\description{
\code{sparsity()} estimates the percentage of entries in a given data matrix \code{D}
whose values are "practically zero". If the absolute value of an entry is
below a given threshold parameter \code{thresh}, then that value is determined
to be "practically zero", increasing the estimated sparsity of \code{D}.
}
\examples{
sparsity(matrix(rep(c(1, 0), 8), 4, 4))
sparsity(matrix(0:8, 3, 3))
sparsity(matrix(0, 3, 3))
}
\seealso{
\code{\link[=matrix_rank]{matrix_rank()}}
}
