\name{ezLev}
\Rdversion{1.1}
\alias{ezLev}
\title{
Function to re-order the levels of a factor
}
\description{
This function provides easy re-ordering of factors.
}
\usage{
ezLev(
	x
	, new_order
)
}
\arguments{
  \item{x}{
Factor vector to be re-ordered.
}
  \item{new_order}{
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

#reorder the levels of flanker
levels(ANT$flanker)
ANT$flanker = ezLev(
	x = ANT$flanker
	, new_order = c('Incongruent Flanker','Neutral Flanker','Congruent Flanker')
)
levels(ANT$flanker)

#alternatively, you can use the level number
levels(ANT$cue)
ANT$cue = ezLev(
	x = ANT$cue
	, new_order = c(4,3,2,1)
)
levels(ANT$cue)


}
