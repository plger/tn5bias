#' getFragmentEndBias
#'
#' Computes the Tn5 bias at the insertion site of each end (start and end) of
#' each fragment. The bias is computed one chromosome at a time for all of its
#' positions (see \code{\link{getSequenceTn5Bias}}, which is fast and
#' multi-threaded) and then looked up at the fragment ends, which is much
#' faster than computing the bias around each end separately unless there are
#' very few fragments.
#'
#' The fragment ends are taken to be the insertion sites, i.e. the fragments
#' are expected to be already shifted (as is the case of the fragment files
#' produced by e.g. 10x Genomics' cellranger-atac), and the bias is the shifted
#' bias (`shift=TRUE`) at those positions.
#'
#' @param fragments Either the path to a tabix-indexed fragment file (BED-like,
#'   with 0-based starts and end-exclusive ends, as in the fragment files), or a
#'   \code{\link[GenomicRanges]{GRanges}} in which the start and end are
#'   themselves the (1-based) positions of the two insertion sites.
#' @param genome A \link[Rsamtools]{FaFile} or \code{BSgenome} from which to
#'   fetch the sequences. Sequence names are matched flexibly to those of the
#'   fragments (e.g. "chr1" and "1").
#' @param bias An optional 10-mer bias table (see \code{\link{getTn5Bias}}).
#' @param scale The bias scale (see \code{\link{estimateBiasScale}}).
#' @param log2 Logical; whether to return the bias as log2(odds ratio) rather
#'   than odds ratio (default).
#' @param chromosomes Optional names of the chromosomes (as in `fragments`) to
#'   restrict the computation to.
#' @param nthreads The number of threads used to compute the bias.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return For a fragment file, a data.frame with one row per fragment, with
#'   the columns `seqnames`, `start`, `end` (as in the file, i.e. 0-based
#'   start), `barcode` (if present in the file), `bias_start` and `bias_end`.
#'   Fragments are in the order of the file within each chromosome, and
#'   chromosomes in the order of the file's index. For a \code{GRanges}, the
#'   same object (same order) with the columns `bias_start` and `bias_end`
#'   added. Ends falling outside the sequence have `NA` bias. Internally the
#'   bias is stored as 100*log2(OR) integers, so values are rounded to
#'   0.01 on the log2 scale (<0.7% on the odds ratio).
#' @export
#' @importFrom Rsamtools seqnamesTabix scanTabix
#' @importFrom GenomeInfoDb seqlengths
#' @examples
#' set.seed(1)
#' s <- paste(sample(c("A","C","G","T"), 1000, replace=TRUE), collapse="")
#' fa <- tempfile(fileext=".fa"); writeLines(c(">chrT", s), fa)
#' Rsamtools::indexFa(fa)
#' fr <- GRanges("chrT", IRanges(c(100, 300), c(180, 450)))
#' getFragmentEndBias(fr, Rsamtools::FaFile(fa))
getFragmentEndBias <- function(fragments, genome, bias=NULL, scale=1,
                               log2=FALSE, chromosomes=NULL, nthreads=1L,
                               verbose=TRUE){
  stopifnot(is(genome, "FaFile") || is(genome, "BSgenome"))
  sl <- seqlengths(genome)
  isGR <- is(fragments, "GRanges")
  if(isGR){
    chrs <- unique(as.character(seqnames(fragments)))
  }else{
    stopifnot(is.character(fragments), length(fragments)==1)
    chrs <- seqnamesTabix(fragments)
  }
  if(!is.null(chromosomes)) chrs <- intersect(chrs, chromosomes)
  gmap <- .matchSeqNames(chrs, names(sl))
  if(anyNA(gmap)){
    if(verbose) message("Skipping sequences not found in the genome: ",
                        paste(chrs[is.na(gmap)], collapse=", "))
    chrs <- chrs[!is.na(gmap)]
  }
  if(length(chrs)==0) stop("No sequence in common between the fragments and ",
                           "the genome.")

  # bias of the 100*log2(OR) integers -> requested scale
  conv <- function(v){
    v <- v/100
    if(log2) v else 2^v
  }

  # positions of the two ends (1-based), for a set of fragments
  res <- lapply(setNames(chrs, chrs), function(chr){
    gn <- gmap[[chr]]
    L <- as.integer(sl[[gn]])
    if(verbose) message("Processing ", chr, "...")
    if(isGR){
      i <- which(as.character(seqnames(fragments))==chr)
      p1 <- start(fragments)[i]; p2 <- end(fragments)[i]
      info <- i
    }else{
      lines <- unlist(scanTabix(fragments, param=GRanges(chr, IRanges(1L,
                                              max(L, 1L)))), use.names=FALSE)
      if(length(lines)==0) return(NULL)
      d <- scan(text=lines, what=list(NULL, 0L, 0L, ""), sep="\t",
                quote="", flush=TRUE, fill=TRUE, quiet=TRUE)
      info <- data.frame(seqnames=chr, start=d[[2]], end=d[[3]],
                         stringsAsFactors=FALSE)
      if(any(nzchar(d[[4]]))) info$barcode <- d[[4]]
      p1 <- d[[2]] + 1L; p2 <- d[[3]]
      rm(lines, d)
    }
    gr <- GRanges(gn, IRanges(1L, L))
    sq <- if(is(genome, "FaFile")) scanFa(genome, gr) else
      Biostrings::getSeq(genome, gr)
    names(sq) <- chr
    b <- getSequenceTn5Bias(sq, bias=bias, nthreads=nthreads, asInteger=TRUE,
                            shift=TRUE, scale=scale, verbose=FALSE)[[1]]
    rm(sq)
    b <- as.integer(b)
    look <- function(p){
      v <- rep(NA_integer_, length(p))
      ok <- !is.na(p) & p >= 1L & p <= L
      v[ok] <- b[p[ok]]
      conv(v)
    }
    out <- list(info=info, bias_start=look(p1), bias_end=look(p2))
    rm(b); gc(FALSE)
    out
  })
  res <- res[!vapply(res, is.null, logical(1))]

  if(isGR){
    fragments$bias_start <- NA_real_; fragments$bias_end <- NA_real_
    for(r in res){
      fragments$bias_start[r$info] <- r$bias_start
      fragments$bias_end[r$info] <- r$bias_end
    }
    return(fragments)
  }
  out <- do.call(rbind, lapply(res, function(r){
    cbind(r$info, bias_start=r$bias_start, bias_end=r$bias_end)
  }))
  row.names(out) <- NULL
  out
}
