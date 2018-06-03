#' Linear interpolation method to calculate y value at time t0.
#'
#' This function uses linear interpolation method to calculate corresponding y value at time t0.
#' 
#' @param x vector of ages.
#' @param y vector of measurements.
#' @param id factor of subject identifiers.
#' @param t0 target age value.
#' @details For each individual's measurements in growth data set, the first time measurement is usually not the same,
#'          people are interested in some measurement index of certain age, for example height of age 6 years old 
#'          for each id. This function can help us compute corresponding measurement at interested age.
#' @return a data frame including id and corresponding y value at t0.
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, Man-Yu Wong \email{mamywong@@ust.hk}
#' @examples
#' require(sitar)
#' data(heights)
#' x <- heights$age
#' h <- heights$height
#' id <- heights$id
#' ## height at age 9
#' h9 <- lim(x,h,id,9)  
#' @export lim
lim <- function(x, y, id, t0){
  linear_interpolation <- function(x0, x, y){
    if(x0 == 0) {
      y0 = 0
      return (y0)
    }
    else{
      if(min(x) >= x0){
        im <- order(x) 
        xs <- x[im]
        ys <- y[im]
        y0 <- (xs[2] - x0) / (xs[2] - xs[1]) * ys[1] - (xs[1] - x0) / (xs[2] - xs[1]) * ys[2]
      }
      else if(max(x) <= x0){
        im <- order(x)
        xs <- x[im]
        ys <- y[im]
        ns <- length(xs)
        x1 <- xs[ns - 1]
        x2 <- xs[ns]
        y1 <- ys[ns - 1]
        y2 <- ys[ns]
        y0 <- (x0 - x1) / (x2 - x1) * (y2 - y1) + y1
      }
      else{
        diff <- x-x0
        ind0 <- 1 * (diff<0)
        ind0 <- which(ind0 == 0)
        ind2 <- ind0[1]
        ind1 <- ind2-1;
        y0 <- (x0 - x[ind1]) / (x[ind2] - x[ind1]) * (y[ind2] - y[ind1]) + y[ind1] 
      }
      return (y0)
    }
  }
  unid <- matrix(unique(id), ncol=1)
  y0m <- apply(unid, 1, function(x1){
    index <- id == x1
    x.id <- x[index]
    y.id <- y[index]
    yx0.id <- linear_interpolation(t0, x.id, y.id)
    lid <- sum(index)
    resul <- c(lid, yx0.id)
    resul
  })
  yx0 <- rep(y0m[2,], y0m[1,])
  newdata <- data.frame(id, yx0)
  newdata
}