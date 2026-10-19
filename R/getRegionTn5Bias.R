#' getRegionTn5Bias
#'
#' Computes an aggregate Tn5 bias score for each region, based on the 
#' per-nucleotide (shifted, i.e. insertion-site) bias of 
#' \code{\link{getSequenceTn5Bias}}. 
#'
#' With default settings (see the `method` argument), the region's overall
#' bias is the mean bias over its positions.
#'
#' @param regions A \code{\link[GenomicRanges]{GRanges}} of regions.
#' @param genome A \link[Rsamtools]{FaFile} or \code{BSgenome} from which to 
#'   fetch the sequences. Sequence names must match those of `regions`.
#' @param scale The bias scale (see \code{\link{estimateBiasScale}}).
#' @param method Either "mean" (default; arithmetic mean of the odds ratios), 
#'   "geometric", "max" (the highest per-position bias in the region) or 
#'   "quantile" (a trimmed maximum: the `q` quantile of the per-position bias).
#' @param q The quantile(s) used when `method="quantile"` (default 0.98).
#' @param log2 Logical; whether to return the score as log2 rather than as 
#'   (observed/expected) odds ratio.
#' @param maxNFraction The maximum fraction of unknown bases (N, or lying 
#'   outside of the sequence) in a region for it to be scored (default 0.1); 
#'   the score of other regions is `NA`. Within scored regions, positions with 
#'   an unknown base within 9bp (i.e. in either of the two 10-mers on which 
#'   their bias is based) are ignored, since their bias is not informative.
#' @param bias An optional 10-mer bias table (see \code{\link{getTn5Bias}}).
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A numeric vector of the same length as `regions` (or, if several 
#'   methods and/or quantiles are requested, a matrix with one column per 
#'   score, e.g. "mean", "max", "quantile_0.98"). `NA` is returned for regions 
#'   with too many unknown bases (see `maxNFraction`), including those 
#'   entirely outside the sequences or on unknown sequences.
#' @export
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom BiocGenerics start end width
#' @importFrom Rsamtools scanFa
#' @examples
#' library(GenomicRanges)
#' set.seed(1)
#' s <- paste(sample(c("A","C","G","T"), 2000, replace=TRUE), collapse="")
#' f <- tempfile(fileext=".fa")
#' Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(chr1=s)), f)
#' Rsamtools::indexFa(f)
#' regions <- GRanges("chr1", IRanges(c(100,500,1000), width=200))
#' getRegionTn5Bias(regions, Rsamtools::FaFile(f))
getRegionTn5Bias <- function(regions, genome, scale=1, 
                             method=c("mean","geometric","max","quantile"), log2=FALSE, 
                             q=0.98, maxNFraction=0.1, bias=NULL, verbose=TRUE){
  if(missing(method)) method <- "mean"
  method <- match.arg(method, several.ok=TRUE)
  stopifnot(is.numeric(q), length(q)>=1, q>=0, q<=1, 
            is.numeric(maxNFraction), length(maxNFraction)==1)
  # one score per column
  cols <- unlist(lapply(unique(method), function(m){
    if(m=="quantile") paste0("quantile_", q) else m
  }))
  stopifnot(is(regions, "GRanges"))
  pad <- 20L # flank needed for the bias to be exact at the edges
  sl <- seqlengths(genome)
  sn <- as.character(seqnames(regions))
  rs <- start(regions); re <- end(regions)
  # part of the region within the sequence
  L <- unname(sl[sn])
  ok <- !is.na(L) & re >= 1L & rs <= L & width(regions) > 0
  if(verbose && !all(ok)) message(sum(!ok), " regions are outside the genome.")
  out <- matrix(NA_real_, nrow=length(regions), ncol=length(cols), 
                dimnames=list(NULL, cols))
  fin <- function(out){
    if(!log2) out <- 2^out
    if(ncol(out)==1) out <- unname(out[,1])
    out
  }
  if(!any(ok)) return(fin(out))
  rs2 <- pmax(rs[ok], 1L); re2 <- pmin(re[ok], L[ok])
  xs <- pmax(rs2-pad, 1L); xe <- pmin(re2+pad, L[ok])
  ext <- GRanges(sn[ok], IRanges(xs, xe))
  sq <- if(is(genome, "FaFile")) scanFa(genome, ext) else 
    Biostrings::getSeq(genome, ext)
  b <- getSequenceTn5Bias(sq, bias=bias, scale=scale, shift=TRUE, 
                          nthreads=1L, verbose=FALSE)
  off <- rs2 - xs + 1L
  len <- re2 - rs2 + 1L
  # The bias of a position is only informative if the sequence is known around
  # it: a position is valid if there is no unknown base (N, or beyond the 
  # sequence's end) within 9 bp of it, i.e. in either of the two 10-mers.
  # Regions with more than `maxNFraction` of unknown bases are not scored.
  nFrac <- rep(NA_real_, length(b))
  out[ok,] <- t(vapply(seq_along(b), function(i){
    known <- charToRaw(as.character(sq[[i]])) %in% charToRaw("ACGTacgt")
    n <- length(known)
    # padded so that missing flanks count as unknown
    cu <- c(0L, cumsum(!c(rep(FALSE, 9L), known, rep(FALSE, 9L))))
    nunk <- cu[seq_len(n)+19L] - cu[seq_len(n)]  # unknown within +/-9 bp
    idx <- seq(off[i], length.out=len[i])
    nFrac[i] <<- (sum(!known[idx]) + (width(regions)[ok][i] - len[i])) /
      width(regions)[ok][i]
    v <- b[[i]][idx][nunk[idx]==0L]
    if(nFrac[i] > maxNFraction || length(v)==0L) return(rep(NA_real_, length(cols)))
    unlist(lapply(unique(method), function(m){
      switch(m, mean=log2(mean(v)), geometric=mean(log2(v)), max=log2(max(v)),
             quantile=log2(quantile(v, q, names=FALSE)))
    }))
  }, numeric(length(cols))))
  if(verbose && any(is.na(out[ok,1]))) 
    message(sum(is.na(out[ok,1])), " regions were not scored (too many unknown bases).")
  fin(out)
}
