\name{ezBoot}
\Rdversion{1.1}
\alias{ezBoot}
\title{
Function to compute bootstrap resampled predictions for each cell in a specified experimental design.
}
\description{
This function is used to compute bootstrap resampled predictions for each cell in a specified experimental design, using either cell means or mixed effects modelling to obtain predictions. The results can be visualized using \code{\link{ezPlotBoot}}.}
\usage{
ezBoot(
    data
    , dv
    , wid
    , within = NULL
    , between = NULL
    , iterations = 1e3
    , lmer = TRUE
    , family = 'gaussian'
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
  \item{iterations}{
Numeric value specifying the number of bootstrap iterations to complete.
}
  \item{lmer}{
Boolean specifying whether to obtain predictions via mixed effects modelling (TRUE; default) or via cell means (FALSE).
}
  \item{family}{
When obtaining predictions via mixed effects modelling (i.e. when \code{lmer=TRUE}), you must specify the residuals family. While the bootstrap is in theory non-parametric, it may be more powerful if you specify a family that might reasonably expected to match your data. For example, if the data are binary outcomes (eg. accuracy), then use the "binomial" family. See \code{\link[lme4]{lmer}}.
}
}
\details{
While \code{within} and \code{between} are both optional, at least one column of \code{data} must be provided to either \code{within} or \code{between}. Any numeric or character variables in \code{data} that are specified as either \code{wid}, \code{within} or \code{between} will be converted to a factor with a warning. Prior to running, \code{dv} is collapsed to a mean for each cell defined by the combination of \code{wid}, \code{within} or \code{between}.
}
\value{
A list containing either two or three components:
\item{fit}{If predictions are obtained by mixed effects modelling, an \code{link[lme4]{lmer}} object consisting of the original mixed effects model}
\item{cells}{A data frame containing predictions for each cell of the design.}
\item{boots}{A data frame containing predictions for each cell of the design from each iteration of the bootstrap procedure.}
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\seealso{
	\code{\link[lme4]{lmer}}, \code{\link{ANT}}, \code{\link{ANT2}}, \code{\link{ezANOVA}}, \code{\link{ezBoot}}, \code{\link{ezCor}}, \code{\link{ezDesign}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}, \code{\link{ezPlotBoot}}, \code{\link{ezPrecis}}, \code{\link{ezStats}}
}
\examples{
#Read in the ANT data (see ?ANT).
data(ANT)
head(ANT)
ezPrecis(ANT)

#Run ezBoot on the accurate RT data
er = ezBoot(
    data = ANT
    , dv = .(error)
    , wid = .(subnum)
    , within = .(cue,flank)
    , between = .(group)
    , family = 'binomial'
    , iterations = 1e1 #1e3 or higher is best for publication
)

#convert the error rate predictions from the logit scale (the default binomial scale) to the % scale
er$cells$value = 100*plogis(er$cells$value)
er$boots$value = 100*plogis(er$boots$value)

#plot the full design
from_ezPlotBoot = ezPlotBoot(
	from_ezBoot = er
	, x = .(flank)
	, split = .(cue)
	, col = .(group)
)
print(from_ezPlotBoot$plot)

#plot the effect of group across the flank*cue design
from_ezPlotBoot = ezPlotBoot(
	from_ezBoot = er
	, x = .(flank)
	, split = .(cue)
	, diff = .(group)
)
print(from_ezPlotBoot$plot)

}
