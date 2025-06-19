#' getTn5Bias
#' 
#' Returns a table of Tn5 bias per 10-mers, obtained from SELMA (Hu et al. 2022)
#'  using ATAC-seq on naked human DNA.
#'
#' @param convert Logical; whether to convert into an (observed/expected) odds 
#'  ratio (default TRUE). If FALSE, the bias values as stored, i.e. rounded, 
#'  100*log2(OR) to ease the memory footprint.
#' @param includeNs Logical; whether to include kmers with leading or trailing 
#'  Ns (default TRUE).
#'
#' @return A data.frame with the columns kmer and bias
#' @export
#'
#' @references 
#' Hu, S.S., Liu, L., Li, Q. et al. Intrinsic bias estimation for improved 
#' analysis of bulk and single-cell chromatin accessibility profiles using 
#' SELMA. Nat Commun 13, 5533 (2022). https://doi.org/10.1038/s41467-022-33194-z
#' 
#' @examples
#' bias <- getTn5Bias()
#' head(bias)
getTn5Bias <- function(convert=TRUE, includeNs=TRUE){
  data("tn5bias")
  kmers <- .getKmers(includeNs=includeNs)
  d <- data.frame(kmer=kmers, bias=head(tn5bias, length(kmers)))
  if(convert) d$bias <- 2^(d$bias/100)
  d
}

.getKmers <- function(includeNs=TRUE){
  nts <- c("A","C","G","T")
  kmers <- expand.grid(nts,nts,nts,nts,nts,nts,nts,nts,nts,nts)
  kmers <- do.call(paste0, as.list(as.data.frame(kmers)))
  if(!includeNs) return(kmers)
  kmers2 <- lapply(1:5, \(i){
    k2 <- do.call(expand.grid, lapply((i+1):10, \(x) nts))
    paste0(paste0(rep("N",i), collapse=""), 
           do.call(paste0, as.list(as.data.frame(k2))))
  })
  kmers3 <- lapply(6:10, \(i){
    k2 <- do.call(expand.grid, lapply(1:(i-1L), \(x) nts))
    paste0(do.call(paste0, as.list(as.data.frame(k2))),
           paste0(rep("N",i-5L), collapse=""))
  })
  c(kmers, unlist(kmers2), unlist(kmers3))
}