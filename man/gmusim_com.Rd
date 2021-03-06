% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gmusim_com.R
\name{gmusim_com}
\alias{gmusim_com}
\title{Fit a growth curve model using the Shape Invariant Model with random effects involved in 
time independent covariates having combined effects.}
\usage{
gmusim_com(x, y, z, w, p1, p2, n, id, df, knots, len, cal = FALSE,
  fixed = random, random = "a+b+c", bounds = 0.04, start, bstart,
  verbose = FALSE, correlation = NULL, weights = NULL, subset = NULL,
  method = "ML", na.action = na.fail, control = nlmeControl(returnObject =
  TRUE))
}
\arguments{
\item{x}{vector of ages.}

\item{y}{vector of measurements.}

\item{z}{data.frame of covariates only affect subject's size.}

\item{w}{data.frame of covariates only affect subject's growth rate.}

\item{p1}{number of columns in data.frame z.}

\item{p2}{number of columns in data.frame w.}

\item{n}{length of x.}

\item{id}{factor of subject identifiers.}

\item{df}{degrees of freedom for cubic regression spline.}

\item{knots}{vector of values for knots (default df quantiles of x distribution).}

\item{len}{time (age) measures after extension for computing aphv, if it is missing, then
it equals round(365*diff(range(x))); when cal=FALSE, it is useless.}

\item{cal}{control for whether to calculate pv, apv and height at pv or not (default FALSE).}

\item{fixed}{character string specifying a, b, c fixed effects (default random).}

\item{random}{character string specifying a, b, c random effects (default "a+b+c").}

\item{bounds}{span of x for regression spline, or fractional extension of range (default 0.04).}

\item{start}{optional numeric vector of initial estimates for the fixed effects, or list of initial
estimates for the fixed and random effects (see \code{\link[nlme]{nlme}}).}

\item{bstart}{optional starting value for fixed effect b, if it is missing, then it equals to mean(x) when calculating.}

\item{verbose}{optional logical value to print information on the evolution of the iterative algorithm (see \code{\link[nlme]{nlme}}).}

\item{correlation}{optional corStruct object describing the within-group correlation structure (see \code{\link[nlme]{nlme}}).}

\item{weights}{optional varFunc object or one-sided formula describing the within-group heteroscedasticity structure (see \code{\link[nlme]{nlme}}).}

\item{subset}{optional expression indicating the subset of the rows of data that should be used
in the fit (see \code{\link[nlme]{nlme}}).}

\item{method}{character string, either "REML" or "ML" (default) (see \code{\link[nlme]{nlme}}).}

\item{na.action}{function for when the data contain NAs (see \code{\link[nlme]{nlme}}).}

\item{control}{list of control values for the estimation algorithm (see \code{\link[nlme]{nlme}}).}
}
\value{
An object inheriting from class gmusim representing the nonlinear mixed-effects model fit, with
all the components returned by nlme (see \code{\link[nlme]{nlmeObject}} for a full description) plus the following
components:

 \strong{bstart:} the value of arg bstart.

 \strong{call.gmusim: } the internal gmusim_com call that produced the object.

 \strong{fitnlme: }  the function returning the predicted value of y.

 \strong{ns: } the lm object providing starting values for the B-spline curve.

 \strong{calindex: } APHV(age at phv), PHV(peak height velocity) and HPHV(height at phv). 

 \strong{fitted.values: } a data frame, including expanded x, fitted values (y) for corresponding expanded x and id. 
 
Generic functions such as print, plot, anova and summary can be used to show the results of the
fit. The functions resid, coef, fitted, fixed.effects, random.effects, predict, getData,
getGroups, getCovariate and getVarCov can be used to extract some of gmusim_com's components.
}
\description{
The basis of the Shape Invariant Model is that a population has a common 
characteristic curve or function, which by shifting and scaling can be made to   
have the form of any individual curve. This function can deal 
with time independent covariates having combined effects, that is, some covariates only have size effects, 
some only have growth rate effects and some others have both size and growth rate effects.
}
\details{
\strong{Start} is an initial estimation vector for fixed effect parameters, it is suggested that the 
initial values for a, b, c are 0, mean(x) and 0, respectively. Note that a, b, c are corresponding to $alpha_0$,
$beta_0$ and -$beta_1$ in model (7) of Beath (2007). And if a covariate has both size and growth rate effect, this
covariate should be in both z and w data.frame. One method for improving the initial 
guess for coefficients of fixed effects is to fit the model without random effects to the pooled data
using a standard nonlinear least squares package. \strong{bstart} allows the origin
of b to be varied. Changing the original of b affects its random effect variance.
}
\examples{
require(sitar)
data(heights)
x <- heights$age
y <- heights$height
id <- heights$id
men <- heights$men
df <- 5
#since negative value are censored, here use absolute value
z <- data.frame(z1=abs(men))
#generate assumed bmi at 6 years old for these 12 girls
set.seed(100)
bmi6 <- runif(12,min=16,max=26)
numid <- as.numeric(table(id))
bmi6 <- rep(bmi6,numid)
w <- data.frame(w1=bmi6)
p1 <- p2 <- 1
n <- length(x)
####Do not calculate phv, aphv and hphv
resu1 <- gmusim_com(x,y,z,w,p1,p2,n,id,df)
summary(resu1)
####Calculate phv, aphv and height at phv (hphv)
resu2 <- gmusim_com(x,y,z,w,p1,p2,n,id,df,cal=TRUE)
summary(resu2)
aphv <- resu2$calindex
fitted.values <- resu2$fitted.values
}
\references{
Beath KJ. Infant growth modelling using a shape invariant model with random effects.
Statistics in Medicine 2007;26:2547-2564.

Cole TJ, Donaldson MD, Ben-Shlomo Y. SITAR--a useful instrument for growth
curve analysis. Int J Epidemiol 2010;39:1558-66.
}
\seealso{
\code{\link[sitar]{sitar}},\code{\link{gmusim}},\code{\link{gmusim_both}},\code{\link{gmusim_gr}}, \code{\link{gmusim_size}}.
}
\author{
Zhiqiang Cao \email{zcaoae@connect.ust.hk}, Man-Yu Wong \email{mamywong@ust.hk}
}
\keyword{models}
\keyword{nonlinear}
\keyword{package}
\keyword{regression}
