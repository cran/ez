\name{ezStats}
\Rdversion{1.1}
\alias{ezStats}
\title{
Function to obtain descriptive statistics from a factorial experiment
}
\description{
This function provides easy computation of descriptive statistics (between-Ss means, between-Ss SD, Fisher's Least Significant Difference) for data from factorial experiments, including purely within-Ss designs (a.k.a. "repeated measures"), purely between-Ss designs, and mixed within-and-between-Ss designs. 
}
\usage{
ezStats(
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
	, between_full = NULL
	, collapse_within = FALSE
)
}
\arguments{
  \item{data}{
Data frame containing the data to be analyzed.
}
  \item{dv}{
.() object specifying the column in \code{data} that contains the dependent variable. Values in this column should be of the numeric class.
}
  \item{sid}{
.() object specifying the column in \code{data} that contains the variable specifying the case/Ss identifier. Values in this column will be converted to factor class if necessary.
}
  \item{within}{
Optional .() object specifying the column(s) in \code{data} that contains independent variables that are manipulated within-Ss. Values in this column will be converted to factor class if necessary.
}
  \item{between}{
Optional .() object specifying the column(s) in \code{data} that contains independent variables that are manipulated between-Ss. Values in this column will be converted to factor class if necessary.
}
  \item{between_full}{
Same as \code{between}, but must specify the full set of between-Ss variables if \code{between} specifies only a subset of the design.
}
  \item{collapse_within}{
Optional boolean to trigger the collapse of a single 2-level within-Ss varbiable to a difference score (useful when obtaining statistics for mixed within-and-between-Ss effects).
}
}
\details{
While \code{within} and \code{between} are both optional, at least one column of \code{data} must be provided to either \code{within} or \code{between}. Any numeric or character variables in \code{data} that are specified as either \code{sid}, \code{within} or \code{between} will be converted to a factor with a warning. Fisher's Least Significant Difference is computed as sqrt(2)*qt(.975,DFd)*sqrt(MSd/N), where N is taken as the mean N per group in cases of unbalanced designs.
}
\value{
A data frame containing the descriptive statistics for the requested effect. N = number of Ss per cell. Mean = between-Ss mean. SD = between-Ss SD. FLSD = Fisher's Least Significant Difference.}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\section{Warning }{
The descriptives include Fisher's Least Significant Difference for the requested effect. In the context of purely within-Ss or purely between-Ss this value may be used for post-hoc multiple comparisons. Note however that in the context of mixed within-and-between-Ss designs, this value can only be used for within-Ss comparisons.
}

\seealso{
\code{\link{ezANOVA}}, \code{\link{ezCor}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}
}
\examples{
#Read in the ANT data (see ?ANT).
data(ANT)

#Show summaries of the ANT data.
head(ANT)
str(ANT)
summary(ANT)

#Compute some useful statistics per cell.
cell_stats = ddply(
	.data = ANT
	, .variables = .( sid , group , cue , flanker )
	, .fun <- function(x){
		#Compute error rate as percent.
		error_rate = (1-mean(x$acc))*100
		#Compute mean RT (only accurate trials).
		mean_rt = mean(x$rt[x$acc==1])
		#Compute SD RT (only accurate trials).
		sd_rt = sd(x$rt[x$acc==1])
		return(c(error_rate=error_rate,mean_rt=mean_rt,sd_rt=sd_rt))
	}
)

#Run an ANOVA on the mean_rt data.
mean_rt_anova = ezANOVA(
	data = cell_stats
	, dv = .(mean_rt)
	, sid = .(sid)
	, within = .(cue,flanker)
	, between = .(group)
)

#Show the ANOVA & assumption tests.
print(mean_rt_anova)

#Compute descriptives for the main effect of group.
group_descriptives = ezStats(
	data = cell_stats
	, dv = .(mean_rt)
	, sid = .(sid)
	, between = .(group)
)

#Show the descriptives.
print(group_descriptives)
}
