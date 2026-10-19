#' getSequenceTn5Bias
#'
#' Computes the Tn5 bias at each position of the input sequence(s), following 
#' the SELMA method: the bias of a cut is the geometric mean of the bias of the 
#' two 10-mers around the two cleavage sites of the Tn5 dimer (9 bp apart).
#'
#' @param x A (set of) DNA sequence(s), either as a \link[Rsamtools]{FaFile} 
#'   object, a character vector, or a DNAString or DNAStringSet.
#' @param bias An optional 10-mer bias table, with the columns 'kmer' and 
#'   'bias'. By default, an (internal) table derived from SELMA will be used 
#'   (see \link{getTn5Bias}).
#' @param nthreads The number of threads to use (default 1). Alternatively, a
#'   \code{BiocParallelParam} object, in which case parallelization is done
#'   across sequences.
#' @param complement Logical; whether to compute the bias of cuts on the 
#'  negative strand (i.e. the 5' end of negative-strand reads, reported at 
#'  their genomic position). Only relevant when `shift=FALSE`.
#' @param asInteger Logical; whether to convert the (observed/expected) odds 
#'  ratio to rounded 100*log2(OR) to ease the memory footprint. Recommended when
#'  computing genome-wide bias.
#' @param shift Logical; whether to shift the bias, so that it represents not the 5'
#'  ends of (positives) fragments, but the insertion site (midpoint between the
#'  cleavages sites on each strand, i.e. 4 bp after the 5' end of a 
#'  positive-strand read). The shifted bias is consistent with the 
#'  shifting often performed on ATAC-seq data, and is thus enabled by default.
#'  Note, however, that this is different from the original SELMA methods.
#'  With `shift=FALSE`, the result is the same as SELMA's.
#' @param scale A numeric factor by which to amplify the (log-scale) bias, i.e.
#'  the bias used is `bias^scale` (default 1). This can be estimated for a given
#'  sample using \code{\link{estimateBiasScale}}.
#' @param verbose Logical; whether to output progress messages.
#'
#' @return An AtomicList of bias per position.
#' @export
#'
#' @importFrom methods is as
#' @importFrom stats setNames
#' @importFrom Biostrings reverseComplement DNAStringSet DNAString width
#' @importFrom Rsamtools scanFa
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam bpnworkers
#' @import IRanges
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom GenomicRanges GRanges
#' @importFrom RcppParallel setThreadOptions RcppParallelLibs
#' @importFrom Rcpp sourceCpp
#' @useDynLib tn5bias, .registration = TRUE
#' @examples
#' fakeSequence <- "NNNNNACGTACGTGCGGGCTATGTCACGTNNNNNNN"
#' getSequenceTn5Bias(fakeSequence)
getSequenceTn5Bias <- function(x, bias=NULL, nthreads=1L, complement=FALSE, 
                               asInteger=FALSE, shift=TRUE, scale=1, verbose=TRUE){
  stopifnot(is.numeric(scale), length(scale)==1, is.finite(scale))
  stopifnot(inherits(nthreads, "BiocParallelParam") || 
              (is.numeric(nthreads) && length(nthreads)==1 && nthreads>=1))
  fa <- NULL
  if(is(x, "DNAString")) x <- DNAStringSet(x)
  if(is(x, "DNAStringSet")){
    totalSize <- sum(width(x))
  }else if(is.character(x)){
    totalSize <- sum(nchar(x))
  }else if(is(x, "FaFile")){
    fa <- x
    sl <- seqlengths(fa)
    totalSize <- sum(sl)
    x <- as.list(sl)
  }else{
    stop("Unrecognized input format!")
  }
  if(is.null(names(x))) names(x) <- paste0("seq", seq_along(x))
  
  if(is.null(bias)){
    bias <- getTn5Bias(convert=TRUE)
  }else{
    stopifnot(c("kmer","bias") %in% colnames(bias))
  }
  # the C++ code indexes kmers alphabetically
  bias <- bias[order(bias$kmer), ]
  bias_vec <- bias$bias^scale
  stopifnot(length(bias_vec) == 4^10)
  
  # Pre-compute RC map once
  all_kmers <- bias$kmer
  rc_kmers <- as.character(reverseComplement(DNAStringSet(all_kmers)))
  rc_map <- match(rc_kmers, all_kmers) - 1L 
  
  if(complement && shift){
    if(verbose) message("The use of shift=TRUE and complement=TRUE together ",
                        "does not make sense, as the two strands are identical",
                        " after shifting. `complement` is set to FALSE.")
    complement <- FALSE
  }
  
  if(length(x)==1) verbose <- FALSE
  if(parOverSeqs <- (length(x) > (totalSize/length(x))/10)){
    if(inherits(nthreads, "BiocParallelParam")){
      bp <- nthreads
    }else if(nthreads>1){
      bp <- MulticoreParam(nthreads, progressbar=verbose)
    }else{
      bp <- SerialParam(progressbar=verbose)
    }
    setThreadOptions(1)
  }else{
    if(inherits(nthreads, "BiocParallelParam")) nthreads <- bpnworkers(nthreads)
    setThreadOptions(nthreads)
    bp <- SerialParam(progressbar=verbose)
  }
  
  # iterate by index, as sequence names are not necessarily unique
  results <- bplapply(seq_along(x), BPPARAM=bp, function(i){
    if(is.null(fa)){
      s_char <- as.character(x[[i]])
    }else{
      name <- names(x)[i]
      s_char <- as.character(scanFa(fa, GRanges(name, IRanges(1L, sl[[name]]))))
    }
    len <- nchar(s_char)
    
    # padding to handle edges (the widest window is 18 bp)
    padded_s <- paste0(strrep("N", 20), s_char, strrep("N", 20))
    seqbias <- get_seq_tn5_bias_par(padded_s, bias_vec, rc_map, shift, 
                                    complement)
    seqbias <- seqbias[21:(20 + len)]
    
    if(asInteger) seqbias <- as.integer(round(100 * log2(seqbias)))
    seqbias
  })
  names(results) <- names(x)
  
  as(results, ifelse(asInteger, "IntegerList", "NumericList"))
}
