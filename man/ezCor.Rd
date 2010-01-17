\name{ezCor}
\Rdversion{1.1}
\alias{ezCor}
\title{
Function to plot a correlation matrix with scatterplots, linear fits, and univariate density plots
}
\description{
This function provides simultaneous visualization of a correlation matrix, scatter-plot with linear fits, and univariate density plots for multiple variables.
}
\usage{
ezCor(
	data
	, r_size_lims = c(10,30)
	, point_alpha = .5
	, density_height = 1
	, density_adjust = 1
	, density_colour = 'white'
	, label_size = 10
	, label_colour = 'black'
	, label_alpha = .5
	, lm_colour = 'red'
	, ci_colour = 'green'
	, ci_alpha = .5	
)
}
\arguments{
  \item{data}{
	Data frame containing named columns of data only.
}
  \item{r_size_lims}{
	Minimum and maximum size of the text reporting the correlation coefficients. Minimum is mapped to coefficients of 0 and maximum is mapped to coefficients of 1, with the mapping proportional to r^2.
}
  \item{point_alpha}{
	Transparency of the data points (1 = opaque).
}
  \item{density_adjust}{
	Adjusts the bandwidth of the univariate density estimator. See "adjust" parameter in \code{\link[stats]{density}}
}
  \item{density_height}{
	Proportion of the facet height taken up by the density plots.
}
  \item{density_colour}{
	Colour of the density plot.
}
  \item{label_size}{
	Size of the variable labels on the diagonal.
}
  \item{label_colour}{
	Colour of the variable labels on the diagonal.
}
  \item{label_alpha}{
	Transparency of the variable labels on the diagonal (1 = opaque).
}
  \item{lm_colour}{
	Colour of the fitted line.
}
  \item{ci_colour}{
	Colour of the confidence interval surrounding the fitted line.
}
  \item{ci_alpha}{
	Transparency of the confidence interval surrounding the fitted line (1 = opaque).
}
}
\value{
A ggplot2 object.
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\seealso{
\code{\link{ezANOVA}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}, \code{\link{ezStats}}
}
\examples{
########
# Set up some fake data
########
library(MASS)
N=100

#first pair of variables
variance1=1
variance2=2
mean1=10
mean2=20
rho = .8
Sigma=matrix(c(variance1,sqrt(variance1*variance2)*rho,sqrt(variance1*variance2)*rho,variance2),2,2)
pair1=mvrnorm(N,c(mean1,mean2),Sigma,empirical=TRUE)

#second pair of variables
variance1=10
variance2=20
mean1=100
mean2=200
rho = -.4
Sigma=matrix(c(variance1,sqrt(variance1*variance2)*rho,sqrt(variance1*variance2)*rho,variance2),2,2)
pair2=mvrnorm(N,c(mean1,mean2),Sigma,empirical=TRUE)

my_data=data.frame(cbind(pair1,pair2))

########
# Now plot
########
p = ezCor(
	data = my_data
)
print(p)

#you can modify the default colours of the correlation coefficients as follows
p = p + scale_colour_manual(values = c('red','blue'))
print(p)
#see the following for alternatives:
# http://had.co.nz/ggplot2/scale_manual.html
# http://had.co.nz/ggplot2/scale_hue.html
# http://had.co.nz/ggplot2/scale_brewer.html

}
