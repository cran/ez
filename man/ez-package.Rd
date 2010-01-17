\name{ez-package}
\Rdversion{1.1}
\alias{ez}
\docType{package}
\title{
Easy analysis and visualization of factorial experiments.
}
\description{
This package facilitates easy analysis of factorial experiments, including purely within-Ss designs (a.k.a. "repeated measures"), purely between-Ss designs, and mixed within-and-between-Ss designs. The functions in this package provide easy access to descriptive statistics (\code{\link{ezStats}}), ANOVA (\code{\link{ezANOVA}}), non-parametric permutation tests (\code{\link{ezPerm}}), and visualization of results (\code{\link{ezPlot}}).}
\details{
\tabular{ll}{
Package: \tab ez\cr
Type: \tab Package\cr
Version: \tab 1.6\cr
Date: \tab 2010-01-17\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
This package contains several useful functions:
\itemize{
\item{\code{\link{ezANOVA}}}{Provides simple interface to ANOVA, including assumption checks.}
\item{\code{\link{ezCor}}}{Function to plot a correlation matrix with scatterplots, linear fits, and univariate density plots}
\item{\code{\link{ezPerm}}}{Provides simple interface to the Permutation test.}
\item{\code{\link{ezPlot}}}{Uses the \code{ggplot2} graphing package to generate plots for any given user-requested effect, by default producing error bars that facilitate visual post-hoc multiple comparisons.}
\item{\code{\link{ezStats}}}{Provides between-Ss descriptive statistics for any given user-requested effect.}
}
}
\author{
Author & Maintainer: Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\seealso{
\code{\link{ezANOVA}}, \code{\link{ezCor}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}, \code{\link{ezStats}}
}
\keyword{ package }
