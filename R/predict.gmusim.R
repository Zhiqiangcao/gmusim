#' Predict the Shape Invariant Model with random effects involved in covariates.
#'
#' Predict method for gmusim objects, based on predict.lme.
#' 
#' @param object an object inheriting from class gmusim.
#' @param newdata an optional data frame to be used for obtaining the predictions. It requires named columns for x, z  
#'        and for id if level = 1, matching the names in object. Note that values of covariates z should be centralizated. 
#'        By default their values are set to the mean so when level = 0 the prediction represents the mean curve. 
#'        Note that factors are coded as instrumental variables for each level, corresponding to the fixed effect 
#'        coefficients in the model, so their names need the level appending..
#' @param level an optional integer giving the level of grouping to be used in obtaining the predictions, 
#'        level 0 corresponding to the population predictions. Defaults to level 1.
#' @param ... other optional arguments, including na.action and naPattern.
#' @details Note that if level = 1, this function calculates predicton for every measurment of individuals; if
#'          level = 0, it calculates mean value of individuals' measurements.  
#' @return A vector of the predictions.
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, Man-Yu Wong \email{mamywong@@ust.hk}
#' @importFrom nlme getData ranef getGroups nlmeControl
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
#' ## fit sitar model with covariates
#' resu1 <- gmusim_both(x,y,z,p,n,id,df)
#' ## predictions at level = 0
#' on.exit(detach(resu1))
#' eval(parse(text = "attach(resu1)"))
#' predict(resu1, newdata=data.frame(x=6:15,z1=rep(mean(men),10)), level=0)
#' ## predictions at level = 1 for all subjects
#' newd <- data.frame(x=heights$age,z1=abs(men)-mean(abs(men)),id=id)
#' on.exit(detach(resu1))
#' eval(parse(text = "attach(resu1)"))
#' fitted.values <- predict(resu1, newdata=newd, level=1)
#' @export 
predict.gmusim <- function(object, newdata, level = 1, ...){
  if (missing(newdata)) 
    newdata <- getData(object)
  oc <- object$call.gmusim
  if (!is.null(newdata$x)) 
    x <- newdata$x
  else newdata$x <- x <- eval(oc$x, newdata)
  abcset <- FALSE
  abc <- ranef(object)
  if (abcset || level == 0) 
    newdata$id <- id <- factor(rep.int(getGroups(object)[1], 
                                       nrow(newdata)))
  else if (!is.null(newdata$id)) 
    id <- factor(newdata$id)
  else newdata$id <- id <- factor(eval(oc$id, newdata))
  on.exit(detach(object))
  eval(parse(text = "attach(object)"))
  fitnlme <- object$fitnlme
  if (is.null(fitnlme)) 
    stop("could not find function \"fitnlme\": please update model")
  argnames <- names(formals(fitnlme))
  args <- setNames(vector("integer", length = length(argnames)), 
                   argnames)
  args <- args[!match(argnames, names(newdata), 0)]
  newdata <- data.frame(newdata, t(args))
  class(object) <- class(object)[-1]
  pred <- predict(object = object, newdata = newdata, 
                       level = level, ...)
  attributes(pred) <- NULL
  return(pred)
}