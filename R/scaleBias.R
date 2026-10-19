#' scaleBias
#'
#' Applies a bias scale (see \code{\link{estimateBiasScale}}) to a 
#' pre-computed bias, which is much faster than re-computing it from the 
#' sequence with the `scale` argument of \code{\link{getSequenceTn5Bias}}.
#'
#' @param x The bias, as returned by \code{\link{getSequenceTn5Bias}} computed 
#'   with `scale=1`, i.e. a \code{NumericList} of odds ratios, an 
#'   \code{IntegerList} of 100*log2(odds ratio) (`asInteger=TRUE`), or a 
#'   numeric/integer vector or list of the same.
#' @param scale The scale, e.g. as estimated by \code{\link{estimateBiasScale}}.
#'   The scaled bias is `x^scale` (or, for integer-encoded values, `x*scale`, 
#'   rounded).
#'
#' @return The scaled bias, of the same type as `x`.
#' @export
#' @examples
#' b <- getSequenceTn5Bias("ACGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGATCGAT", 
#'                         verbose=FALSE)
#' scaleBias(b, 1.5)
#' # equivalent to (but much faster than):
#' getSequenceTn5Bias("ACGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGATCGAT", scale=1.5,
#'                    verbose=FALSE)
scaleBias <- function(x, scale){
  stopifnot(is.numeric(scale), length(scale)==1, is.finite(scale))
  if(is.list(x) || is(x, "List")){
    isInt <- if(is(x, "List")) is(x, "IntegerList") else 
      all(vapply(x, is.integer, logical(1)))
    if(is(x, "List")) return(.scaleBiasList(x, scale, isInt))
    return(lapply(x, scaleBias, scale=scale))
  }
  if(is.integer(x)) return(as.integer(round(x*scale)))
  x^scale
}

.scaleBiasList <- function(x, scale, isInt){
  if(isInt) as(round(as(x, "NumericList")*scale), "IntegerList") else x^scale
}
