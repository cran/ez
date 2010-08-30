\name{ezANOVA}
\Rdversion{1.1}
\alias{ezANOVA}
\title{
Function to perform a factorial ANOVA
}
\description{
This function provides easy analysis of data from factorial experiments, including purely within-Ss designs (a.k.a. "repeated measures"), purely between-Ss designs, and mixed within-and-between-Ss designs, yielding ANOVA results and assumption checks.
}
\usage{
ezANOVA(
    data
    , dv
    , wid
    , within = NULL
    , between = NULL
    , observed = NULL
    , diff = NULL
    , reverse_diff = FALSE
)
}
\arguments{
  \item{data}{
Data frame containing the data to be analyzed.
}
  \item{dv}{
.() object specifying the column in \code{data} that contains the dependent variable. Values in this column must be numeric.
}
  \item{wid}{
.() object specifying the column in \code{data} that contains the variable specifying the case/Ss identifier.
}
  \item{within}{
Optional .() object specifying one or more columns in \code{data} that contain predictor variables that are manipulated (or observed) within-Ss.
}
  \item{between}{
Optional .() object specifying one or more columns in \code{data} that contain predictor variables that are manipulated (or observed) between-Ss.
}
  \item{observed}{
Optional .() object specifying one or more columns in \code{data} that contain predictor variables that are observed variables (i.e. not manipulated). The presence of observed variables affects the computation of the generalized eta-squared measure of effect size reported by \code{\link{ezANOVA}}.
}
  \item{diff}{
Optional .() object specifying a 2-level within-Ss varbiable to collapse to a difference score.
}
  \item{reverse_diff}{
Optional boolean to trigger reversal of the difference collapse requested by \code{diff}.
}
}
\details{
While \code{within} and \code{between} are both optional, at least one column of \code{data} must be provided to either \code{within} or \code{between}. Any numeric or character variables in \code{data} that are specified as either \code{wid}, \code{within} or \code{between} will be converted to a factor with a warning. Prior to running, \code{dv} is collapsed to a mean for each cell defined by the combination of \code{wid}, \code{within} or \code{between}.
}
\value{
A list containing one or more of the following components:
\item{ANOVA}{A data frame containing the ANOVA results.}
\item{Mauchly's Test for Sphericity}{If any within-Ss variables with >2 levels are present, a data frame containing the results of Mauchly's test for Sphericity. Only reported for effects >2 levels because sphericity necessarily holds for effects with only 2 levels.}
\item{Sphericity Corrections}{If any within-Ss variables are present, a data frame containing the Greenhouse-Geisser & Huynh-Feldt epsilon values, and corresponding corrected p-values.}
\item{Levene's Test for Homgeneity}{If the design is purely between-Ss, a data frame containing the results of Levene's test for Homgeneity of variance.}
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}

\section{Warning}{
The statistical computing in \code{ezAnova()} is driven by the \code{\link[car]{Anova}} function from the \code{car} package, using the univariate Type-II test. If there are too few Ss in the data set for \code{\link[car]{Anova}} to compute its MANOVA, \code{ezAnova()} will revert to using \code{\link[stats]{aov}} for computing the ANOVA, in which case no assumption tests are provided. When assumption tests are provided via \code{\link[car]{Anova}}, Huynh-Feldt corrected p-values where the Huynh-Feldt epsilon >1 will use 1 as the correction epsilon.
}
\note{
 Some column names in the output data frames are abbreviated to conserve space:
\itemize{
\item{DFn}{Degrees of Freedom in the numerator (a.k.a. DFeffect).}
\item{DFd}{Degrees of Freedom in the denominator (a.k.a. DFerror).}
\item{SSn}{Sum of Squares in the numerator (a.k.a. SSeffect).}
\item{SSd}{Sum of Squares in the denominator (a.k.a. SSerror).}
\item{F}{F-value.}
\item{p}{p-value (probability of the null hypothesis given the data).}
\item{p<.05}{Highlights p-values less than the traditional alpha level of .05.}
\item{ges}{Generalized Eta-Squared (see in references below: Bakeman, 2005) equal to SSn divided by: The effect's SSn plus all unique SSd's in the ANOVA plus all SSn's for effects involving observed variables minus SSn for the effect if the effect involves an observed variable. }
\item{GGe}{Greenhouse-Geisser epsilon.}
\item{p[GGe]}{p-value after correction using Greenhouse-Geisser epsilon.}
\item{p[GGe]<.05}{Highlights p-values (after correction using Greenhouse-Geisser epsilon) less than the traditional alpha level of .05.}
\item{HFe}{Huynh-Feldt epsilon.}
\item{p[HFe]}{p-value after correction using Huynh-Feldt epsilon.}
\item{p[HFe]<.05}{Highlights p-values (after correction using Huynh-Feldt epsilon) less than the traditional alpha level of .05.}
\item{W}{Mauchly's W statistic}
}
}
\references{
Bakeman, R. (2005). Recommended effect size statistics for repeated measures designs. Behavior Research Methods, 37 (3), 379-384.
}
\seealso{
	\code{\link{ANT}}, \code{\link{ANT2}}, \code{\link{ezANOVA}}, \code{\link{ezBoot}}, \code{\link{ezCor}}, \code{\link{ezDesign}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}, \code{\link{ezPlotBoot}}, \code{\link{ezPrecis}}, \code{\link{ezStats}}
}
\examples{
#Read in the ANT data (see ?ANT).
data(ANT)
head(ANT)
ezPrecis(ANT)


#Run an ANOVA on the mean correct RT data.
mean_rt_anova = ezANOVA(
    data = ANT[ANT$error==0,]
    , dv = .(rt)
    , wid = .(subnum)
    , within = .(cue,flank)
    , between = .(group)
)

#Show the ANOVA & assumption tests.
print(mean_rt_anova)

#Run an ANOVA on the mean_rt data, ignoring group.
mean_rt_anova2 = ezANOVA(
    data = ANT[ANT$error==0,]
    , dv = .(rt)
    , wid = .(subnum)
    , within = .(cue,flank)
)

#Show the ANOVA & assumption tests.
print(mean_rt_anova2)

#Run a purely between-Ss ANOVA on the mean_rt data.
mean_rt_anova3 = ezANOVA(
    data = ANT[ANT$error==0,]
    , dv = .(rt)
    , wid = .(subnum)
    , between = .(group)
)

#Show the ANOVA & assumption tests.
print(mean_rt_anova3)
}
