\name{ezMixed}
\Rdversion{1.1}
\alias{ezMixed}
\title{
Provides tests for each fixed effect in an mixed effects model
}
\description{
This function permits specification of a mixed effects model with one or more simple random effects and one or more multiple fixed effects, yielding tests of each fixed effect and interactions thereof. Mixed models are built using the \code{\link[lme4]{lmer}} function and compared to their respective null model via the \code{\link[stats]{anova}} function.
}
\usage{
ezMixed(
	data
	, dv
	, random
	, fixed
	, family = gaussian
	, alarm = TRUE
	, depth = 0
	, return_models = FALSE
	, return_anovas = FALSE
)
}
\arguments{
  \item{data}{
Data frame containing the data to be analyzed.
}
  \item{dv}{
.() object specifying the column in \code{data} that contains the dependent variable. Values in this column must be numeric.
}
  \item{random}{
.() object specifying one or more columns in \code{data} that contain random effects.
}
  \item{fixed}{
.() object specifying one or more columns in \code{data} that contain fixed effects.
}
  \item{family}{
 a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. See \code{\link[stats]{family}}.
}
  \item{alarm}{
Logical. If TRUE (default), call the \code{\link[utils]{alarm}} function when \code{\link{ezMixed}} completes.
}
  \item{depth}{
Integer specifying the depth of interaction between fixed effects to test. The default value (0) tests the full depth, but can be time-intensive for complex models.
}
  \item{return_models}{
Logical. If TRUE, the returned list object will also include each lmer model (can become memory intensive for complex models).
}
  \item{return_anovas}{
Logical. If TRUE, the returned list object will also include the anova object produced by each test.
}
}
\details{
Each test necessitates building two models, a model that contains the effect of interest plus lower order effects and a model that contains only the lower order effects. For example, a test of the effect of a variable V1 will compare the models \code{dv~(1|random)+V1} versus \code{dv~(1|random)+V1}, a test of the interaction between variables V1 and V2 will compare the model \code{dv~(1|random)+V1+V2+V1:V2} versus \code{dv~(1|random)+V1+V2}, etc. Building each model can be computationally intensive (especially when the family isn't gaussian), so results of each test are printed to screen while the function progresses.
}
\value{
A list whose first element is a summary of the nested comparisons, including p-values from the likelihood ratio test and the likelihood ratio. Depending on the \code{return_models} and \code{return_anovas} arguments, the returned list may also include the fitted models for each effect and the anova object on which the summary is based.
}
\author{
Michael A. Lawrence \email{Mike.Lawrence@dal.ca}
}
\seealso{
    \code{\link[lme4]{lmer}}, \code{\link{ANT}}, \code{\link{ANT2}}, \code{\link{ezANOVA}}, \code{\link{ezBoot}}, \code{\link{ezBootPlot}}, \code{\link{ezCor}}, \code{\link{ezDesign}}, \code{\link{ezMixed}}, \code{\link{ezPerm}}, \code{\link{ezPlot}}, \code{\link{ezPrecis}}, \code{\link{ezPredict}}, \code{\link{ezStats}}
}
\examples{
#Read in the ANT data (see ?ANT).
data(ANT)
head(ANT)
ezPrecis(ANT)

#Run ezMixed on the accurate RT data
rt = ezMixed(
    data = ANT[ANT$error==0,]
    , dv = .(rt)
    , random = .(subnum)
    , fixed = .(cue,flank,group)
)
print(rt$summary)

#Run ezMixed on the error rate data
er = ezMixed(
    data = ANT
    , dv = .(error)
    , random = .(subnum)
    , fixed = .(cue,flank,group)
    , family = 'binomial'
)
print(er$summary)

}
