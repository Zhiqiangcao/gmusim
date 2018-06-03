#' Do interpolation for original time measurements and covariates. 
#'
#' This function evenly partitions time measurements and produces additional ones for an 
#' individual in his/her age range. It then uses \code{\link{predict.gmusim}} to obtain additional fitted 
#' values to calculate age at peak velocity (apv), peak velocity (pv) and height (or weight) at peak velocity (ypv).
#'
#' @param x vector of ages (assume x is centered).
#' @param z data.frame of time independent covariates (assume z is centered).
#' @param p number of columns in data.frame z.
#' @param id factor of subject identifiers.
#' @param idmat matrix of unique id, note that this matrix has been setted to 1 column. 
#' @param n time (age) measures after extension for computing aphv (default is round(365*diff(range(x)))).
#' @details For some individuals, the number of measurements is small. In order to calculate accurate apv 
#'          (age at peak velocity), pv (peak velocity) and ypv (height at peak velocity or weight at peak 
#'          velocity), it is necessary to perform interpolation for original time measurements and obtain 
#'          additional predictions for each individual. This function extends time (age) measurements 
#'          between min(time) and max(time), as well as covariates z, note that output of z has been centralized.  
#' @return a data frame including expanded x(age), centralization covariates z and corresponding id.
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, Man-Yu Wong \email{mamywong@@ust.hk}
#' @examples
#' require(sitar)
#' data(heights)
#' x <- heights$age
#' h <- heights$height
#' id <- heights$id
#' men <- heights$men
#' men <- abs(men)-mean(abs(men))
#' z <- data.frame(z1=men)
#' p <- 1
#' idmat <- matrix(unique(id), ncol = 1)
#' newdata <- exdataz(x, z, p, id, idmat)
#' @export 
exdataz <- function(x, z, p, id, idmat, n=round(365 * diff(range(x)))){
  exdatainn <- function(x, z, p, id, muid, n){
    zz <- paste("z", 1:p, sep = "")
    npt <- n/diff(range(x))
    extage <- apply(idmat, 1, function(x1){
      index <- id==x1
      id.x <- x[index]
      nt <- floor(npt*diff(range(id.x))) + 1
      id.z <- z[index,]
      newx <- seq(min(id.x), max(id.x), length = nt)
      if(p == 1){
        newz <- matrix(rep(id.z[1], nt), ncol = p)
      }
      else {
        newz <- matrix(rep(id.z[1,], nt), byrow = T, ncol = p)
      }
      newid <- rep(x1,nt)
      extx <- data.frame(x = newx, z = newz, id = newid)
      colnames(extx) <- c("x", zz, "id")
      extx
    })
    df <- extage[[1]][FALSE, ]
    for(dft in extage) df <- rbind(df, dft)
    df
  }
  newdata <- exdatainn(x, z, p, id, idmat, n)
  ncol <- dim(newdata)[2]
  for(i in 1:ncol){
    if(class(newdata[,i]) == "list") newdata[,i] <- as.numeric(newdata[,i])
  }
  newdata
}