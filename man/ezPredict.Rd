\name{ezPredict}
\Rdversion{1.1}
\alias{ezPredict}
\title{
Computes predicted values from a mixed effects model
}
\description{
This function computes the predicted values from the fixed effects design of a mixed effects model.
}
\usage{
ezPredict( 
	fit
	, to_predict = NULL
	, numeric_res = 1e3
)
}
\arguments{
  \item{fit}{
Fitted \code{\link[lme4]{lmer}} object.
}
  \item{to_predict}{
Optional data frame containing the fixed effects design to predict. If absent, the function will assume that the full design from the provided fitted model is requested.
}
  \item{numeric_res}{
Integer value specifying the sampling resolution of any numeric fixed effect. Has no effect if non-NULL value supplied to \code{to_predict}. If \code{to_predict} is null and a numeric fixed effect is encountered in the fitted model, then predictions will be obtained at this many evenly spaced intervals between the minimum and maximum values in the original fitted data.
}
}
\value{
A data frame containing the prediction value (and estimated variance of this value) for each cell in the fixed effects design.
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\seealso{
    \code{\link{ANT}}, \code{\link{ANT2}}, \code{\link{ezANOVA}}, \code{\link{ezBoot}}, \code{\link{ezBootPlot}}, \code{\link{ezCor}}, \code{\link{ezDesign}}, \code{\link{ezMixed}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}, \code{\link{ezPredict}}, \code{\link{ezStats}}
}
\examples{
#Read in the ANT data (see ?ANT).
data(ANT)
head(ANT)

#fit a mixed effects model to the error rate data
er_fit = lmer(
	formula = error ~ cue*flank*group + (1|subnum)
	, family = binomial
	, data = ANT
)

#obtain the predictions from the model
er_preds = ezPredict(
	fit = er_fit
)

#compute 95% CI for each prediction
er_preds$lo = er_preds$value - qnorm(.975)*sqrt(er_preds$var)
er_preds$hi = er_preds$value + qnorm(.975)*sqrt(er_preds$var)


#visualize the predictions
ggplot(
	data = er_preds
	, mapping = aes(
		x = flank
		, y = value
		, ymin = lo
		, ymax = hi
	)
)+
geom_point(
	alpha = .75
)+
geom_line(
	alpha = .5
)+
geom_errorbar(
	alpha = .5
)+
facet_grid(
	cue ~ group
)+
labs(
	y = 'Error Rate (log odds)'
)

#collapse incongruent and congruent to a difference score
## Note that the variance of a difference score is the sum 
## of the variances contributing to the score.
er_preds2 = ddply(
	.data = er_preds
	, .variables = .(cue,group)
	, .fun = function(x){
		to_return = data.frame(
			value = x$value[x$flank=='Incongruent']-x$value[x$flank=='Congruent']
			, var = x$var[x$flank=='Incongruent']+x$var[x$flank=='Congruent']
		)
	}
)

#compute the 95% CI for the Incongruent-Congruent difference
er_preds2$lo = er_preds2$value - qnorm(.975)*sqrt(er_preds2$var)
er_preds2$hi = er_preds2$value + qnorm(.975)*sqrt(er_preds2$var)

#visualize the predictions
ggplot(
	data = er_preds2
	, mapping = aes(
		x = cue
		, y = value
		, ymin = lo
		, ymax = hi
	)
)+
geom_point(
	alpha = .75
)+
geom_line(
	alpha = .5
)+
geom_errorbar(
	alpha = .5
)+
facet_grid(
	. ~ group
)+
labs(
	y = 'Error Rate Difference (Incongruent - Congruent, log-odds)'
)

}
