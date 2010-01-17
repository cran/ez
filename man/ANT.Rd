\name{ANT}
\Rdversion{1.1}
\alias{ANT}
\docType{data}
\title{
ANT data
}
\description{
\emph{Simulated} data from then Attention Network Test (see reference below), consisting of 2 within-Ss variables ("cue" and "flanker"), 1 between-Ss variable ("group") and 2 dependent variables (response time, "rt", and accuracy, "acc")
}
\usage{data(ANT)}
\format{
  A data frame with 5760 observations on the following 10 variables.
  \describe{
    \item{\code{sid}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10} \code{11} \code{12} \code{13} \code{14} \code{15} \code{16} \code{17} \code{18} \code{19} \code{20}}
    \item{\code{group}}{a factor with levels \code{Control} \code{Treatment}}
    \item{\code{block}}{a numeric vector}
    \item{\code{trial}}{a numeric vector}
    \item{\code{cue}}{a factor with levels \code{No Cue} \code{Center Cue} \code{Double Cue} \code{Spatial Cue}}
    \item{\code{flanker}}{a factor with levels \code{Neutral Flanker} \code{Congruent Flanker} \code{Incongruent Flanker}}
    \item{\code{location}}{a factor with levels \code{down} \code{up}}
    \item{\code{direction}}{a factor with levels \code{left} \code{right}}
    \item{\code{rt}}{a numeric vector}
    \item{\code{acc}}{a numeric vector}
  }
}
\references{
J Fan, BD McCandliss, T Sommer, A Raz, MI Posner (2002). Testing the efficiency and independence of attentional networks. \emph{Journal of Cognitive Neuroscience}, \bold{14}, 340-347.
}
\seealso{
\code{\link{ezANOVA}}, \code{\link{ezCor}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}, \code{\link{ezStats}}
}
\examples{
data(ANT)
head(ANT)
str(ANT)
summary(ANT)
}
\keyword{datasets}
