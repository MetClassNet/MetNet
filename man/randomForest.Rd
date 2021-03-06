% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistical.R
\name{randomForest}
\alias{randomForest}
\title{Create an adjacency matrix based on random forest}
\usage{
randomForest(x, ...)
}
\arguments{
\item{x}{matrix, where columns are the samples and the rows are features
(metabolites), cell entries are intensity values}

\item{...}{parameters passed to `GENIE3`}
}
\value{
matrix, matrix with the importance of the links inferred from
random forest algorithm implemented by `GENIE3`
}
\description{
`randomForest` infers an adjacency matrix using
random forest using the `GENIE3` function from the
`GENIE3` package. `randomForest` returns the importance of the link
between features in the form of an adjacency matrix.
}
\details{
For use of the parameters used in the `GENIE3` function,
refer to `?GENIE3::GENIE3`. The arguments `regulators` and `targets` are
set to `NULL`. Element \eqn{w_{i,j}} (row i, column j) gives the importance
of the link from i to j.
}
\examples{
data("x_test", package = "MetNet")
x <- x_test[1:10, 3:ncol(x_test)]
x <- as.matrix(x)
randomForest(x)

}
\author{
Thomas Naake, \email{thomasnaake@googlemail.com}
}
