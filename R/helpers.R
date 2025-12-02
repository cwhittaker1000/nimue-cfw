#' Convert a vector of parameter values and a vector of the times it changes into a time-series of the parameter values.
#'
#' This essentially replicates the interpolate(_, "constant") function of odin, but returns a time series of values at the requested times t.
#' The implementation for 1D vector x's is not finished and should be checked before usage, alternatively x can be replaced with indexes from x's time dimension.
#'
#' @param t Vector of times we want parameter values for.
#' @param x Vector (list or matrix or array) of parameter values that correspond to each tt_x.
#' @param tt_x Vector of times at which the parameter x changes.
#' @return A vector of values for parameter x that correspond to the times specified by t
#'
#' @export
block_interpolate <- function(t, x, tt_x){
  indexes <- stats::approx(x = tt_x, y=seq_along(tt_x), xout = t, method = "constant", rule = 2)$y
  if(is.vector(x)){
    x[indexes]
  } else {
    if(is.list(x)){
      val <- x[indexes]
    } else if(is.array(x)){
      val <- x[indexes,,]
    } else if(is.matrix(x)){
      val <- x[indexes,]
    }
    if(length(indexes) == 1){
      val[[1]]
    } else {
      val
    }
  }
}
