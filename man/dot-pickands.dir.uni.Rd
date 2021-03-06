% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/liouville_extreme.R
\name{.pickands.dir.uni}
\alias{.pickands.dir.uni}
\title{Pickands dependence function for the copula domain of attraction of Liouville survival copulas}
\usage{
.pickands.dir.uni(t, alpha, rho)
}
\arguments{
\item{t}{pseudo-angle in (0,1)}

\item{alpha}{vector of Dirichlet allocations. Currently must be of length 2}

\item{rho}{index of regular variation parameter}
}
\value{
value of Pickands function for the scaled extremal Dirichlet  model
}
\description{
Pickands dependence function as in \cite{Belzile (2014), Proposition 40 and Example 4}
Returns the Pickands dependence function of the copula domain of attraction
of the survival copula, the scaled Dirichlet extreme value model. Currently only implemented in the bivariate case.
Setting \code{rho=1} yields the same output as the function in the \code{evd} package.
}
