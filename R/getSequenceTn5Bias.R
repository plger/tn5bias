#' getSequenceTn5Bias
#'
#' @param x A (set of) DNA sequence(s), either as a \link[Rsamtools]{FaFile} 
#'   object, a character vector, or a DNAString or DNAStringSet.
#' @param bias An optional 10-mer bias table, with the columns 'kmer' and 
#'   'bias'. By default, an (internal) table derived from SELMA will be used 
#'   (see \link{getTn5Bias}).
#' @param nthreads The number of threads to use (default 1). Alternatively, 
#' @param complement Logical; whether to compute bias for the complement of `x`.
#'  This is useful to get the bias for the negative strand.
#' @param asInteger Logical; whether to convert the (observed/expected) odds 
#'  ratio to rounded 100*log2(OR) to ease the memory footprint. Recommended when
#'  computing genome-wide bias.
#' @param shift Logical; whether to shift the bias, so that it represents not the 5'
#'  ends of (positives) fragments, but the insertion site (midpoint between the
#'  cleavages sites on each strand). The shifted bias is consistent with the 
#'  shifting often performed on ATAC-seq data, and is thus enabled by default.
#'  Note, however, that this is different from the original SELMA methods.
#' @param verbose Logical; whether to output progress messages.
#'
#' @return An AtomicList of bias per position.
#' @export
#'
#' @importFrom Biostrings reverseComplement DNAStringSet DNAString
#' @importFrom Rsamtools scanFa
#' @importFrom BiocParallel bplapply MulticoreParam bpnworkers
#' @import IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom RcppParallel setThreadOptions RcppParallelLibs
#' @importFrom Rcpp sourceCpp
#' @useDynLib tn5bias, .registration = TRUE
#' @examples
#' fakeSequence <- "NNNNNACGTACGTGCGGGCTATGTCACGTNNNNNNN"
#' getSequenceTn5Bias(fakeSequence)
getSequenceTn5Bias <- function(x, bias=NULL, nthreads=1L, complement=FALSE,
                               asInteger=FALSE, shift=TRUE, verbose=TRUE){
  stopifnot(inherits(nthreads, "BiocParallelParam") || 
              (is.numeric(nthreads) && length(nthreads)==1 && nthreads>=1))
  if(is(x, "DNAString")) x <- DNAStringSet(x)
  if(is(x, "DNAStringSet")){
    totalSize <- sum(width(x))
  }else if(is.character(x)){
    totalSize <- sum(nchar(x))
  }else if(is(x, "FaFile")){
    x <- lapply(seqlengths(x), FUN=\(l) list(file=x, length=l))
    totalSize <- sum(seqlengths(x))
  }else{
    stop("Unrecognized input format!")
  }
  if(is.null(names(x))) names(x) <- paste0("seq", seq_along(x))
  
  if(is.null(bias)){
    bias <- getTn5Bias(convert=TRUE)
  }else{
    stopifnot(c("kmer","bias") %in% colnames(bias))
  }
  
  if(complement && shift){
    if(verbose) message("The use of shift=TRUE and complement=TRUE together ",
                        "does not make sense, as the two strands are identical",
                        " after shifting. `complement` is set to FALSE.")
    complement <- TRUE
  }
    
  
  if(length(x)==1) verbose <- FALSE
  if(parOverSeqs <- (length(x) > (totalSize/length(x))/10)){
    if(inherits(nthreads, "BiocParallelParam")){
      bp <- nthreads
    }else if(nthreads>1){
      bp <- BiocParallel::MulticoreParam(nthreads, progress=verbose)
    }else{
      bp <- BiocParallel::SerialParam(progress=verbose)
    }
    setThreadOptions(1)
  }else{
    if(inherits(nthreads, "BiocParallelParam")) nthreads <- bpnworkers(nthreads)
    setThreadOptions(nthreads)
    bp <- BiocParallel::SerialParam(progress=verbose)
  }
  
  as(bplapply(setNames(names(x), names(x)), BPPARAM=bp, \(name){
    s <- x[[name]]
    if(is(s, "DNAString")){
      if(complement) s <- Biostrings::reverseComplement(s)
      s <- as.character(s)
    }else if(is(s, "character")){
      if(complement) s <- as.character(reverseComplement(DNAString(s)))
    }else{
      s <- scanFa(s$file, GRanges(name, IRanges(1L, s$length)))
      if(complement) s <- Biostrings::reverseComplement(s)
      s <- as.character(s)
    }
    size <- nchar(s)
    # if(FALSE && verbose){
    #   prop <- round(100*computedSize/totalSize)
    #   message(prop, "% done... now running sequence ", name)
    # }
    s <- paste0("NNNNN",s,"NNNNN")
    if(parOverSeqs){
      seqbias <- get_seq_tn5_bias(s, bias)
    }else{
      seqbias <- get_seq_tn5_bias_par(s, bias)
    }
    # remove the flanking bases
    seqbias <- seqbias[-c(seq_len(5), seq(from=size+6L, to=size+10L, by=1L))]
    if(shift) seqbias <- c(rep(1,4), head(seqbias, size-4L))
    if(asInteger){
      seqbias <- as.integer(round(100*log2(seqbias)))
    }
    if(complement) seqbias <- rev(seqbias)
    seqbias
  }), ifelse(asInteger, "IntegerList", "NumericList"))
}
