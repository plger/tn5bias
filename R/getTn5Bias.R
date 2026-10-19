#' getTn5Bias
#' 
#' Returns a table of Tn5 bias per 10-mers, obtained from SELMA (Hu et al. 2022)
#'  using ATAC-seq on naked human DNA.
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
  d <- d[order(d$kmer), ]
  row.names(d) <- NULL
  d
}

.getKmers <- function(){
  nts <- c("A", "C", "G", "T")
  kmers <- expand.grid(nts, nts, nts, nts, nts, nts, nts, nts, nts, nts, 
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  do.call(paste0, kmers)
}
