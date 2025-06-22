#' getMotifTn5Bias
#' 
#' Returns the per-position Tn5 bias of a DNA motif. Note that this is very 
#' slow and doesn't scale well to larger motifs.
#'
#' @param x A DNA motif, in any format recognized by 
#'   \link[universalmotif]{convert_motifs}.
#' @param bias An optional 10-mer bias table, with the columns 'kmer' and 
#'   'bias'. By default, an (internal) table derived from SELMA will be used 
#'   (see \link{getTn5Bias}).
#' @param ... Any argument passed to \code{\link{getSequenceTn5Bias}}.
#'
#' @return A vector of bias per position.
#' @export
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom matrixStats rowProds
#' @importFrom universalmotif convert_motifs 
#' @examples
#' # generate a dummy motif example
#' mo <- matrix(c(0,0,0,100,100,0,0,0,0,60,40,0,0,100,0,0,0,0,100,0,0,50,50,0),
#'              nrow=4)
#' row.names(mo) <- c("A","C","G","T")
#' getMotifTn5Bias(mo)
getMotifTn5Bias <- function(x, bias=NULL, ...){
  x <- convert_motifs(x)@motif
  if(is.null(bias)) bias <- getTn5Bias()
  x <- t(t(x)/colSums(x))
  choices <- lapply(seq_len(ncol(x)), function(i){
    factor(row.names(x)[which(x[,i] >= 0.1)], row.names(x))
  })
  kmers <- do.call(expand.grid, choices)
  prob <- vapply(seq_len(ncol(kmers)), FUN = function(i){
    x[as.integer(kmers[, i]),i]
  }, FUN.VALUE = numeric(nrow(kmers)))
  prob <- matrixStats::rowProds(prob)
  b <- getSequenceTn5Bias(DNAStringSet(do.call(paste0, kmers)), ...)
  b <- matrix(unlist(b), nrow=length(b))
  colSums(b*prob)/sum(prob)
}