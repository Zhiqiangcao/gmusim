#' Fit a growth curve model using the Shape Invariant Model with random effects involved in time independent covariates affecting 
#' subject's growth rate. 
#'
#' The basis of the Shape Invariant Model is that a population has a common 
#' characteristic curve or function, which by shifting and scaling can be made to   
#' have the form of any individual curve. This function can deal 
#' with time independent covariates affecting subject's growth rate, note that every covariate 
#' in data.frame z only has growth rate effect.
#'
#' @param x vector of ages.
#' @param y vector of measurements.
#' @param z data.frame of time independent covariates.
#' @param p number of columns in data.frame z.
#' @param n length of x.
#' @param id factor of subject identifiers.
#' @param df degrees of freedom for cubic regression spline.
#' @param knots vector of values for knots (default df quantiles of x distribution).
#' @param len time (age) measures after extension for computing aphv, if it is missing, then
#'        it equals round(365*diff(range(x))); if cal=FALSE, it is useless. 
#' @param cal control for whether to calculate pv, apv and height at pv or not (default FALSE). 
#' @param fixed character string specifying a, b, c fixed effects (default random).
#' @param random character string specifying a, b, c random effects (default "a+b+c").
#' @param bounds span of x for regression spline, or fractional extension of range (default 0.04).
#' @param start optional numeric vector of initial estimates for the fixed effects, or list of initial
#'        estimates for the fixed and random effects (see \code{\link[nlme]{nlme}}).
#' @param bstart optional starting value for fixed effect b, if it is missing, then it equals to mean(x) when calculating.
#' @param verbose optional logical value to print information on the evolution of the iterative algorithm (see \code{\link[nlme]{nlme}}).
#' @param correlation optional corStruct object describing the within-group correlation structure (see \code{\link[nlme]{nlme}}).
#' @param weights optional varFunc object or one-sided formula describing the within-group heteroscedasticity structure (see \code{\link[nlme]{nlme}}).
#' @param subset optional expression indicating the subset of the rows of data that should be used
#'        in the fit (see \code{\link[nlme]{nlme}}).
#' @param method character string, either "REML" or "ML" (default) (see \code{\link[nlme]{nlme}}).
#' @param na.action function for when the data contain NAs (see \code{\link[nlme]{nlme}}).
#' @param control list of control values for the estimation algorithm (see \code{\link[nlme]{nlme}}).
#' @details \strong{Start} is an initial estimation vector for fixed effect parameters, it is suggested that the 
#' initial values for a, b, c are 0, mean(x) and 0, respectively. Note that a, b, c are corresponding to $alpha_0$,
#' $beta_0$ and -$beta_1$ in model (7) of Beath (2007). One method for improving the initial 
#' guess for coefficients of fixed effects is to fit the model without random effects to the pooled data
#' using a standard nonlinear least squares package. \strong{bstart} allows the origin
#' of b to be varied. Changing the original of b affects its random effect variance. 
#' @return An object inheriting from class gmusim representing the nonlinear mixed-effects model fit, with
#' all the components returned by nlme (see \code{\link[nlme]{nlmeObject}} for a full description) plus the following
#' components:
#' 
#'  \strong{bstart:} the value of arg bstart.
#' 
#'  \strong{call.gmusim: } the internal gmusim_gr call that produced the object.
#' 
#'  \strong{fitnlme: }  the function returning the predicted value of y.
#' 
#'  \strong{ns: } the lm object providing starting values for the B-spline curve.
#' 
#'  \strong{calindex: } APHV(age at phv), PHV(peak height velocity) and HPHV(height at phv). 
#' 
#'  \strong{fitted.values: } a data frame, including expanded x, fitted values (y) for corresponding expanded x and id. 
#' 
#' Generic functions such as print, plot, anova and summary canbe used to show the results of the
#' fit. The functions resid, coef, fitted, fixed.effects, random.effects, predict, getData,
#' getGroups, getCovariate and getVarCov can be used to extract some of sitar_gr's components.
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, Man-Yu Wong \email{mamywong@@ust.hk}
#' @references Beath KJ. Infant growth modelling using a shape invariant model with random effects.
#' Statistics in Medicine 2007;26:2547-2564.
#' 
#' Cole TJ, Donaldson MD, Ben-Shlomo Y. SITAR--a useful instrument for growth
#' curve analysis. Int J Epidemiol 2010;39:1558-66.
#' @export
#' @seealso \code{\link[sitar]{sitar}},\code{\link{gmusim}},\code{\link{gmusim_both}},\code{\link{gmusim_com}},\code{\link{gmusim_size}}.
#' @keywords package nonlinear regression models
#' @examples
#' require(sitar)
#' data(heights)
#' x <- heights$age
#' y <- heights$height
#' men <- heights$men
#' id <- heights$id
#' df <- 5
#' #since negative value are censored, here use absolute value
#' z <- data.frame(z1=abs(men))
#' p <- 1
#' n <- length(x)
#' ## Do not calculate phv, aphv and hphv
#' resu1 <- gmusim_gr(x,y,z,p,n,id,df)
#' summary(resu1)
#' ## Calculate phv, aphv and height at phv (hphv)
#' resu2 <- gmusim_gr(x,y,z,p,n,id,df,cal=TRUE)
#' summary(resu2)
#' aphv <- resu2$calindex
#' fitted.values <- resu2$fitted.values
gmusim_gr <- function(x, y, z, p, n, id, df, knots, len, cal=FALSE, 
                    fixed = random, random = "a+b+c", bounds = 0.04, start,
                    bstart, verbose = FALSE, correlation = NULL, weights = NULL, 
                    subset = NULL, method = "ML", na.action = na.fail,
                    control = nlmeControl(returnObject = TRUE)){
  diff.quot <- function(x, y) {
    ## Difference quotient (central differences where available)
    n <- length(x); i1 <- 1:2; i2 <- (n-1):n
    c(diff(y[i1]) / diff(x[i1]), (y[-i1] - y[-i2]) / (x[-i1] - x[-i2]),
      diff(y[i2]) / diff(x[i2]))
  }
  mcall <- match.call()
  if (missing(z))
    stop("one data.frame must be specified")
  if (missing(df) & missing(knots)) 
    stop("either df or knots must be specified")
  if (!missing(df) & !missing(knots)) 
    cat("both df and knots specified - df redefined from knots\n")
  if (missing(knots)) 
    knots <- quantile(x, (1:(df - 1))/df)
  else df <- length(knots) + 1
  if (length(x) <= df) 
    stop("too few data to fit spline curve")
  if (length(bounds) == 1) 
    bounds <- range(x) + abs(bounds) * c(-1, 1) * diff(range(x)) 
  if (length(bounds) != 2) 
    stop("bounds should be length 1 or 2")
  if (missing(bstart)) bstart <- mean(x)
  x <- x- bstart
  knots <- knots - bstart      
  bounds <- bounds - bstart 
  zm <- matrix(rep(apply(z, 2, mean), n), byrow = TRUE, ncol = p)
  z <- z - zm 
  zz <- paste("z", 1:p, sep = "")
  bz <- paste("bz", 1:p, sep = "")
  spline.lm <- lm(y ~ ns(x, knots = knots, Boundary.knots = bounds))
  coefpar <- coef(spline.lm)
  if (missing(start)) start <- c(coefpar[c(2:(df + 1))], rep(0, p), coefpar[1])
  fix <- fixed
  if (!grepl("a", fix)) 
    fix <- paste("a", fix, sep = "+")
  ss <- paste("s", 1:df, sep = "")
  fixed <- c(ss, bz)
  pars <- c("x", zz, ss, bz)
  names(model) <- model <- letters[1:3]
  if (is.null(subset)) subset <- 1:length(x)   
  fulldata <- data.frame(x, y, z, id, subset)
  colnames(fulldata) <- c("x","y", zz, "id", "subset")
  for (l in model){
    if (!grepl(l, fix) && !grepl(l, random)){
      model[l] <- NA
      next
    }
    pars <- c(pars, l)
    if (!grepl(l, fix)) next
    fixed <- c(fixed, l)
    if (l == "b") 
      start <- c(start, 0)
    if (l == "c") 
      start <- c(start, 0)
  }
  pars <- paste(pars, collapse = ",")
  fixed <- paste(fixed, collapse = "+")
  sscomma <- paste(ss, collapse = ",")
  zcomma <- paste(zz, collapse = ",")
  bzcomma <- paste(bz, collapse = ",")
  nsd <- paste(model["a"],"+")
  fitenv <- NULL
  fitcode <- c("fitenv <- new.env()", "fitenv$fitnlme <- function($pars) {", 
               "as.vector( $nsd","(as.matrix(cbind($sscomma)) * as.matrix(ns(",
               "(exp((as.matrix(cbind($bzcomma)) * as.matrix(cbind($zcomma))) %*%", 
               "matrix(rep(1,p), ncol=1)) * x - (b)) * exp(c),", 
               "knots=knots, Boundary.knots=bounds))) %*%", "matrix(rep(1,df), ncol=1))", 
               "}", "on.exit(detach(fitenv))", "attach(fitenv)", 
               "nlme(y ~ fitnlme($pars),", "fixed = $fixed ~ 1,", 
               "random = $random ~ 1 | id,", "data = fulldata,", 
               "start = start, correlation = correlation,", 
               "weights = weights, subset = subset, method = method,", 
               "na.action = na.action, control = control, verbose = verbose)")
  for (i in c("random", "pars", "fixed", "sscomma", "zcomma", "bzcomma", "nsd")) 
    fitcode <- gsub(paste("$", i, sep = ""), get(i), fitcode, fixed = TRUE)
  nlme.out <- eval(parse(text = fitcode))
  nlme.out$fitnlme <- fitenv$fitnlme
  nlme.out$bstart <- bstart
  nlme.out$call.gmusim <- mcall 
  nlme.out$ns <- coefpar
  if (!"gmusim_gr" %in% class(nlme.out)) 
    class(nlme.out) <- c("gmusim_gr", class(nlme.out))
  if (cal == TRUE){
    if (missing(len)) len <- round(365 * diff(range(x))) 
    idmat <- matrix(unique(id), ncol = 1)
    newdata <- exdataz(x, z, p, id, idmat, n=len)
    on.exit(detach(nlme.out))
    eval(parse(text = "attach(nlme.out)"))
    fit.y <- predict(nlme.out, newdata=newdata) 
    fit.x <- newdata$x
    fit.id <- newdata$id
    newd1 <- data.frame(fit.x, fit.y, fit.id)
    calindex <- apply(idmat, 1, function(x1){
      ind <- newd1$fit.id == x1
      x.id <- newd1$fit.x[ind]
      y.id <- newd1$fit.y[ind]
      dydx <- diff.quot(x.id, y.id)
      maxid <- which.max(dydx)
      pv <- dydx[maxid]
      apv <- x.id[maxid] 
      hpv <- y.id[maxid] 
      resu <- c(pv, apv, hpv)
    })
    calindex <- cbind(as.numeric(idmat), t(calindex))
    colnames(calindex) <- c("id", "pv", "apv", "hpv")
    calindex <- as.data.frame(calindex)
    calindex$apv <- calindex$apv + bstart
    nlme.out$calindex <- calindex
    nlme.out$fitted.values <- newd1
  }
  nlme.out
}

