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
#' # get an example motif
#' jaspar <- universalmotif::read_jaspar(system.file("extdata", "jaspar.txt",
#'                                                   package="universalmotif"))
#' getMotifTn5Bias(jaspar[[1]])
getMotifTn5Bias <- function(x, bias=NULL, ...){
  x <- convert_motifs(x)@motif
  if(is.null(bias)) bias <- getTn5Bias()
  x <- t(t(x)/colSums(x))
  choices <- lapply(seq_len(ncol(x)), \(i) row.names(x)[which(x[,i]>=0.1)] )
  kmers <- do.call(expand.grid, choices)
  prob <- vapply(seq_len(ncol(kmers)), FUN=\(i) x[as.integer(kmers[,i])],
                 FUN.VALUE=numeric(nrow(kmers)))
  prob <- matrixStats::rowProds(prob)
  b <- getSequenceTn5Bias(DNAStringSet(do.call(paste0, kmers)), ...)
  b <- matrix(unlist(b), nrow=length(b))
  colSums(b*prob)/sum(prob)
}