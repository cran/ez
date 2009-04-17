\name{ez-package}
\Rdversion{1.1}
\alias{ez-package}
\docType{package}
\title{
Easy analysis and visualization of factorial experiments.
}
\description{
This package facilitates easy analysis of factorial experiments, including purely within-Ss designs (a.k.a. "repeated measures"), purely between-Ss designs, and mixed within-and-between-Ss designs. The functions in this package provide easy access to descriptive statistics (\code{ezStats()}), ANOVA (\code{ezANOVA()}, including non-parametric permutation tests), and visualization of results (\code{ezPlot()}).}
\details{
\tabular{ll}{
Package: \tab ez\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-04-26\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
This package contains three useful functions. \code{ezANOVA()} provides simple interface to ANOVA (including assumption checks). \code{ezStats()} provides descriptive statistics for any given user-requested effect. \code{ezPlot()} uses the ggplot2 graphing package to generate plots for any given user-requested effect, by default producing error bars that facilitate visual post-hoc multiple comparisons.
}
\author{
Author & Maintainer: Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\keyword{ package }
