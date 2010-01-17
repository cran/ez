\name{ezPerm}
\Rdversion{1.1}
\alias{ezPerm}
\title{
Function to perform a factorial permutation test
}
\description{
This function provides easy non-parametric permutation test analysis of data from factorial experiments, including purely within-Ss designs (a.k.a. "repeated measures"), purely between-Ss designs, and mixed within-and-between-Ss designs.
}
\usage{
ezPerm(
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
	, perms
)
}
\arguments{
  \item{data}{
Data frame containing the data to be analyzed.
}
  \item{dv}{
.() object specifying the column in \code{data} that contains the dependent variable. Values in this column must be numeric.
}
  \item{sid}{
.() object specifying the column in \code{data} that contains the variable specifying the case/Ss identifier.
}
  \item{within}{
Optional .() object specifying one or more columns in \code{data} that contain independent variables that are manipulated within-Ss.
}
  \item{between}{
Optional .() object specifying one or more columns in \code{data} that contain independent variables that are manipulated between-Ss.
}
  \item{perms}{
An integer \code{>} 0 specifying the number of permutations to compute.
}
}
\details{
While \code{within} and \code{between} are both optional, at least one column of \code{data} must be provided to either \code{within} or \code{between}. Any numeric or character variables in \code{data} that are specified as either \code{sid}, \code{within} or \code{between} will be converted to a factor with a warning. The expected standard deviation of p-values is approximately \code{sqrt(true_p*(1-true_p)/perms)}; significance tests using an alpha of .05 should therefore employ at least 1e3 permutations. As the permutation test is computationally intensive, it is advisable to pre-test smaller values of \code{perms} and extrapolate to estimate the total test duration before attempting a full run. To facilitate such extrapolation, test duration is provided in the output after running a permutation test.
}
\value{
A list containing one or more of the following components:
\item{Permutation Test}{A data frame containing the permutation test results.}
\item{Test Duration}{An estimate of the test duration in seconds.}
}
\section{Warning }{
\code{ezPerm()} is a work in progress. Under the current implementation, only main effects may be trusted.
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\seealso{
\code{\link{ezANOVA}}, \code{\link{ezCor}}, \code{\link{ezPlot}}, \code{\link{ezStats}}
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

#Compute the grand mean RT per Ss.
gmrt = ddply(
	.data = cell_stats
	, .variables = .( sid , group )
	, .fun <- function(x){
		y = mean(x$mean_rt)
		return(c(y=y))
	}
)

#Run a purely between-Ss ANOVA on the mean_rt data. 
# (Completes after ~30s on a 2.4GHz processor).
mean_rt_perm = ezPerm(
	data = gmrt
	, dv = .(y)
	, sid = .(sid)
	, between = .(group)
	, perms = 1e3
)

#Show the Permutation test.
print(mean_rt_perm)
}
