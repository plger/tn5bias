#' getTn5Bias
#' 
#' Returns a table of Tn5 bias per 10-mers, obtained from SELMA (Hu et al. 2022)
#'  using ATAC-seq on naked human DNA.
#'
#' @details Preparing this table takes a fraction of a second. The functions 
#' that need it (e.g. \code{\link{getSequenceTn5Bias}}, 
#' \code{\link{insertionProfiles}}, \code{\link{getRegionTn5Bias}}, 
#' \code{\link{getFragmentEndBias}}) build it internally when their `bias` 
#' argument is NULL. When calling such functions many times (e.g. in a loop or 
#' over many samples), it is therefore faster to call `getTn5Bias()` once and 
#' pass the result as their `bias` argument (`kmerBias` in 
#' \code{\link{estimateBiasScale}}).
#'
#' @param convert Logical; whether to convert into an (observed/expected) odds 
#'  ratio (default TRUE). If FALSE, the bias values as stored, i.e. rounded, 
#'  100*log2(OR) to ease the memory footprint.
#'
#' @return A data.frame with the columns kmer and bias
#' @export
#' @importFrom utils data head
#'
#' @references 
#' Hu, S.S., Liu, L., Li, Q. et al. Intrinsic bias estimation for improved 
#' analysis of bulk and single-cell chromatin accessibility profiles using 
#' SELMA. Nat Commun 13, 5533 (2022). https://doi.org/10.1038/s41467-022-33194-z
#' 
#' @examples
#' bias <- getTn5Bias()
#' head(bias)
getTn5Bias <- function(convert=TRUE){
  tn5bias <- get(data("tn5bias", package="tn5bias", envir=environment()), 
                 envir=environment())
  kmers <- .getKmers()
  d <- data.frame(kmer=kmers, bias=head(tn5bias, length(kmers)))
  if(convert) d$bias <- 2^(d$bias/100)
  d <- d[order(d$kmer, method="radix"), ]
  row.names(d) <- NULL
  d
}

# all 10-mers, in the order of the stored data (the first nucleotide varies 
# fastest)
.getKmers <- function(){
  nts <- c("A", "C", "G", "T")
  kmers <- nts
  for(i in 2:10) kmers <- paste0(rep(kmers, 4), rep(nts, each=length(kmers)))
  kmers
}
