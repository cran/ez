\name{ezPrecis}
\Rdversion{1.1}
\alias{ezPrecis}
\title{
Function to obtain a structure summary of a given data frame.
}
\description{
This function provides a structure summary of a given data frame.
}
\usage{
ezPrecis( 
    data
    , transpose = TRUE
)
}
\arguments{
  \item{data}{
Data frame containing the data to be analyzed.
}
  \item{transpose}{
Logical. If TRUE (default), triggers tranposition of the resulting summary data frame (useful when there are many columns in the original data frame, leading the untransposed summary data frame to wrap).
}
}
\details{
This function was inspired by the whatis() function from the YaleToolkit package.
}
\value{
A data frame containing the descriptive information about each column in the specified data frame:
\item{type}{
This row indicates the type of data R thinks is in each column. Recall that when R imports data to a data frame, each column is given a label that indicates what type of information is in that column (character, numeric, or a factor data).
}
\item{missing}{
This row reports a count of the number of missing values in each column.
}
\item{unique}{
This row reports a count of the number of unique values in each column.
}
\item{min}{
This row reports the minimum value found in each column. If the column data is numeric this is straightforward. If the column data is factored, the first level is reported. If the column data is character, the alphabetically first string is reported.
}
\item{max}{
This row reports the maximum value found in each column. If the column data is numeric this is straightforward. If the column data is factored, the last level is reported. If the column data is character, the alphabetically last string is reported.
}
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\seealso{
    \code{\link{ANT}}, \code{\link{ANT2}}, \code{\link{ezANOVA}}, \code{\link{ezBoot}}, \code{\link{ezBootPlot}}, \code{\link{ezCor}}, \code{\link{ezDesign}}, \code{\link{ezMixed}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}, \code{\link{ezPrecis}}, \code{\link{ezPredict}}, \code{\link{ezStats}}
}
\examples{
#Read in the ANT2 data (see ?ANT2).
data(ANT2)
head(ANT2)

#Show a summary of the ANT2 data.
ezPrecis(ANT2)
}
