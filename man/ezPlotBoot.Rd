\name{ezPlotBoot}
\Rdversion{1.1}
\alias{ezPlotBoot}
\title{
Function to plot bootstrapped predictions and confidence intervals.
}
\description{
This function provides easy visualization of any given user-requested effect from the bootstrap predictions computed by \code{\link{ezBoot}}.
}
\usage{
ezPlotBoot(
    from_ezBoot
    , x
    , split = NULL
    , row = NULL
    , col = NULL
    , do_lines = TRUE
    , bar_width = NULL
    , to_numeric = NULL
    , x_lab = NULL
    , y_lab = NULL
    , split_lab = NULL
    , levels = NULL
    , diff = NULL
    , reverse_diff = FALSE
    , row_y_free = FALSE
)
}
\arguments{
  \item{from_ezBoot}{
An list object resulting from a call to \code{\link{ezBoot}}.
}
  \item{x}{
.() object specifying the variable to plot on the x-axis.
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
  \item{do_lines}{
Logical. If TRUE, lines will be plotted connecting groups of points.
}
  \item{bar_width}{
Optional numeric value specifying custom widths for the error bar hat.
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
  \item{diff}{
Optional .() object specifying a 2-level varbiable to collapse to a difference score.
}
  \item{reverse_diff}{
Optional boolean to trigger reversal of the difference collapse requested by \code{diff}.
}
  \item{row_y_free}{
Logical. If TRUE, then rows will permit different y-axis scales.
}
}
\value{
A list with 4 components:
\item{plot}{A printable/modifiable ggplot2 object.}
\item{cells}{A data frame containing predictions for each cell in the requested design.}
\item{boots}{A data frame containing each iteration of bootstrap predictions, collased to the requested design.}
\item{boot_stats}{A data frame containing the "lo" and "hi" bounds of the 95 percent bootstrap confidence interval of each cell in the requested design.}
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\seealso{
	\code{\link{ANT}}, \code{\link{ANT2}}, \code{\link{ezANOVA}}, \code{\link{ezBoot}}, \code{\link{ezCor}}, \code{\link{ezDesign}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}, \code{\link{ezPlotBoot}}, \code{\link{ezPrecis}}, \code{\link{ezStats}}
}
\examples{
#Read in the ANT data (see ?ANT).
data(ANT)
head(ANT)
ezPrecis(ANT)

#Run ezBoot on the accurate RT data
rt = ezBoot(
    data = ANT
    , dv = .(rt)
    , wid = .(subnum)
    , within = .(cue,flank)
    , between = .(group)
    , iterations = 1e1 #1e3 or higher is best for publication
)

#plot the full design
from_ezPlotBoot = ezPlotBoot(
	from_ezBoot = rt
	, x = .(flank)
	, split = .(cue)
	, col = .(group)
)
print(from_ezPlotBoot$plot)

#plot the effect of group across the flank*cue design
from_ezPlotBoot = ezPlotBoot(
	from_ezBoot = rt
	, x = .(flank)
	, split = .(cue)
	, diff = .(group)
)
print(from_ezPlotBoot$plot)

#plot the flank*cue design, averaging across group
from_ezPlotBoot = ezPlotBoot(
	from_ezBoot = rt
	, x = .(flank)
	, split = .(cue)
)
print(from_ezPlotBoot$plot)

}
