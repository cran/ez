\name{relev}
\Rdversion{1.1}
\alias{relev}
\title{
Function to re-order the levels of a factor
}
\description{
This function provides easy re-ordering of factors.
}
\usage{
relev(
	x
	, levs
)
}
\arguments{
  \item{x}{
Factor vector to be re-ordered.
}
  \item{levs}{
Character vector indicating the new order of levels.
}
}
\details{
See examples.
}
\value{
The re-ordered factor vector.
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\seealso{
\code{\link{relevel}}
}
\examples{
#Read in the ANT data (see ?ANT).
data(ANT)

#show the levels of flanker
levels(ANT$flanker)

#reorder and show the levels of flanker
ANT$flanker = relev(ANT$flanker,c('Incongruent Flanker','Neutral Flanker','Congruent Flanker'))
levels(ANT$flanker)
}
