\name{ezPlot}
\Rdversion{1.1}
\alias{ezPlot}
\title{
Function to plot data from a factorial experiment
}
\description{
This function provides easy visualization of any given user-requested effect from factorial experiments, including purely within-Ss designs (a.k.a. "repeated measures"), purely between-Ss designs, and mixed within-and-between-Ss designs. By default, Fisher's Least Significant Difference is computed to provide error bars that facilitate visual post-hoc multiple comparisons (see Warning section below).
}
\usage{
ezPlot(
	data
	, dv
	, sid
	, within = NULL
	, between = NULL
	, between_full = NULL
	, x
	, do_lines = TRUE
	, do_bars = TRUE
	, bar_width = NULL
	, bar_size = NULL
	, split = NULL
	, row = NULL
	, col = NULL
	, to_numeric = NULL
	, x_lab = NULL
	, y_lab = NULL
	, split_lab = NULL
	, levels = NULL
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
  \item{x}{
.() object specifying the variable to plot on the x-axis.
}
  \item{do_lines}{
Logical. If TRUE, lines will be plotted connecting groups of points.
}
  \item{do_bars}{
Logical. If TRUE, error bars will be plotted.
}
  \item{bar_width}{
Optional numeric value specifying custom widths for the error bar hat.
}
  \item{bar_size}{
Optional numeric value or vector specifying custom size of the error bars.
}
  \item{split}{
Optional .() object specifying a variable by which to split the data into different shapes/colors (and line types, if do_lines==TRUE).
}
  \item{row}{
Optional .() object specifying a variable by which to split the data into rows.
}
  \item{col}{
Optional .() object specifying a variable by which to split the data into columns.
}
  \item{to_numeric}{
Optional .() object specifying any variables that need to be converted to the numeric class before plotting.
}
  \item{x_lab}{
Optional character string specifying the x-axis label.
}
  \item{y_lab}{
Optional character string specifying the y-axis label.
}
  \item{split_lab}{
Optional character string specifying the key label.
}
  \item{levels}{
Optional named list where each item name matches a factored column in \code{data} that needs either reordering of levels, renaming of levels, or both. Each item should be a list containing named elements \code{new_order} or \code{new_names} or both.
}
  \item{collapse_within}{
Optional boolean to trigger the collapse of a single 2-level within-Ss varbiable to a difference score (useful when plotting mixed within-and-between-Ss effects).
}
}
\details{
While \code{within} and \code{between} are both optional, at least one column of \code{data} must be provided to either \code{within} or \code{between}. Any numeric or character variables in \code{data} that are specified as either \code{sid}, \code{within} or \code{between} will be converted to a factor with a warning. Fisher's Least Significant Difference is computed as sqrt(2)*qt(.975,DFd)*sqrt(MSd/N), where N is taken as the mean N per group in cases of unbalanced designs.
}
\value{
A printable/modifiable ggplot2 object.
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}

\section{Warning }{
The default error bars are Fisher's Least Significant Difference for the plotted effect, facilitating visual post-hoc multiple comparisons. Note however that in the context of mixed within-and-between-Ss designs, these bars can only be used for within-Ss comparisons.
}

\seealso{
\code{\link{ezANOVA}}, \code{\link{ezCor}}, \code{\link{ezPerm}}, \code{\link{ezStats}}
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

#Plot the main effect of group.
group_plot = ezPlot(
	data = cell_stats
	, dv = .(mean_rt)
	, sid = .(sid)
	, between = .(group)
	, x = .(group)
	, do_lines = FALSE
	, x_lab = 'Group'
	, y_lab = 'RT (ms)'
)

#Show the plot.
print(group_plot)

#Plot the cue*flanker interaction.
cue_by_flanker_plot = ezPlot(
	data = cell_stats
	, dv = .(mean_rt)
	, sid = .(sid)
	, within = .(cue,flanker)
	, x = .(flanker)
	, split = .(cue)
	, x_lab = 'Flanker'
	, y_lab = 'RT (ms)'
	, split_lab = 'Cue'
)

#Show the plot.
print(cue_by_flanker_plot)

#Plot the group*cue*flanker interaction.
group_by_cue_by_flanker_plot = ezPlot(
	data = cell_stats
	, dv = .(mean_rt)
	, sid = .(sid)
	, within = .(cue,flanker)
	, between = .(group)
	, x = .(flanker)
	, split = .(cue)
	, col = .(group)
	, x_lab = 'Flanker'
	, y_lab = 'RT (ms)'
	, split_lab = 'Cue'
)

#Show the plot.
print(group_by_cue_by_flanker_plot)

#Re-plot the main effect of group, using the levels
##argument to re-arrange/rename levels of group
group_plot = ezPlot(
	data = cell_stats
	, dv = .(mean_rt)
	, sid = .(sid)
	, between = .(group)
	, x = .(group)
	, do_lines = FALSE
	, x_lab = 'Group'
	, y_lab = 'RT (ms)'
	, levels = list(
		group = list(
			new_order = c('Treatment','Control')
			, new_names = c('Treatment\nGroup','Control\nGroup')
		)
	)
)

#Show the plot.
print(group_plot)

}
