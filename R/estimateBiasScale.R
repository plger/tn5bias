#' estimateBiasScale
#'
#' Estimates, for a given sample, the factor `k` by which the (log-scale) 
#' kmer bias of \code{\link{getTn5Bias}} should be amplified (or attenuated) so
#' that it best describes the observed insertions. The bias is assumed to be 
#' the same across samples, but its strength to vary (e.g. with Tn5 
#' concentration or incubation conditions). The corrected bias is then 
#' `bias^k`, which can be obtained from a pre-computed bias with \code{\link{scaleBias}} 
#' (or from the sequence with the `scale` argument of 
#' \code{\link{getSequenceTn5Bias}}).
#' 
#' Windows are tiled along the sequences and, conditional on the total number
#' of insertions in a window, the insertions are modelled as multinomially 
#' distributed across the window's positions with probabilities proportional to 
#' `bias^k`. Because this is conditional on the window total, differences in 
#' accessibility between windows have no influence; only the within-window 
#' distribution matters. This is exact for naked DNA or mitochondrial DNA (whose
#' accessibility is uniform), and a good approximation for windows of a size 
#' within which accessibility is roughly flat. Windows are sub-sampled so as 
#' to be spread evenly across the range of mean bias.
#'
#' @param x Either the path to an indexed BAM file (requires `genome`; the 5' 
#'   ends of all reads are used, so it works for single- as well as paired-end 
#'   data), a tabix-indexed fragment file, a \code{GRanges} of width-1 cut 
#'   sites (see \code{\link{insertionProfiles}}; for these, `seqs` should 
#'   preferably be given), or the insertion counts per position, i.e. a named list (or 
#'   \code{IntegerList}, \code{RleList}) with one vector per sequence. 
#'   These should be the (shifted) insertions, matching the default 
#'   (`shift=TRUE`) output of \code{\link{getSequenceTn5Bias}}, i.e. +4 for 
#'   positive-strand reads and -4 for negative-strand reads, using the 1-based 
#'   positions of read starts and ends respectively (this is done 
#'   automatically for BAM input).
#' @param bias The bias, as returned by \code{\link{getSequenceTn5Bias}} 
#'   (computed with `scale=1`) for the same sequences. Odds ratios are expected 
#'   unless the values are integers (`asInteger=TRUE`). Ignored (computed from 
#'   `genome`) for BAM input.
#' @param seqs Optional names of the sequences to use (e.g. "chrM" or "MT"). 
#'   By default all are used for list input. For BAM input, the mitochondrial 
#'   sequence is used if present (with enough reads), otherwise the sequence 
#'   with the most reads (to limit memory usage, the bias and counts are 
#'   handled for the selected sequences only, so consider specifying a subset 
#'   in the case of large genomes).
#' @param windowSize The window size (default 200).
#' @param minCount The minimum number of insertions for a window to be used.
#' @param nWindows The (maximum) number of windows to use, spread evenly across
#'   the range of mean bias.
#' @param nBins The number of bias bins across which windows are balanced, and 
#'   in which the diagnostics are reported.
#' @param interval The range of scales considered.
#' @param genome A \link[Rsamtools]{FaFile} or \code{BSgenome}; required for BAM 
#'   input, from which the sequences (and thereby the bias) are obtained. 
#'   Sequence names must match those of the BAM file.
#' @param minMapq For BAM input, the minimum mapping quality (default 10).
#' @param nthreads The number of threads used to compute the bias (BAM input).
#' @param kmerBias An optional 10-mer bias table (see \code{\link{getTn5Bias}}),
#'   used to compute the bias from `genome` (BAM/fragment file input).
#'
#' @return A list with the estimated `scale`, its (model-based, which assumes no
#'   overdispersion and is thus optimistic) standard error `se`, the number of 
#'   windows and insertions used, and a `diagnostic` data.frame, with one row 
#'   per bias bin, giving the observed and expected insertions under no scaling
#'   (`expected1`) and under the estimated scale (`expectedK`).
#' @export
#' @importFrom stats uniroot
#' @examples
#' set.seed(1)
#' s <- paste(sample(c("A","C","G","T"), 50000, replace=TRUE), collapse="")
#' b <- getSequenceTn5Bias(s, verbose=FALSE)
#' # simulate insertions from a bias twice as strong as the reference:
#' p <- log2(as.numeric(b[[1]]))
#' y <- rpois(length(p), 0.5 * 2^(2*p) / mean(2^(2*p)))
#' estimateBiasScale(list(seq1=y), b, minCount=10)$scale
estimateBiasScale <- function(x, bias=NULL, seqs=NULL, windowSize=200L, 
                              minCount=20L, nWindows=5000L, nBins=20L,
                              interval=c(-1,6), genome=NULL, minMapq=10L,
                              nthreads=1L, kmerBias=NULL){
  if(is.character(x) || is(x, "GRanges")){
    if(is.character(x) && length(x)!=1) stop("Only one file can be given.")
    if(is.null(genome)) stop("`genome` is required for BAM/fragment file input.")
    xb <- .bamInsertionsAndBias(x, genome, seqs, minCount, minMapq, nthreads,
                                kmerBias)
    x <- xb$x; bias <- xb$bias; seqs <- names(x)
  }
  if(is.null(bias)) stop("`bias` is required.")
  stopifnot(length(x)==length(bias), !is.null(names(x)))
  if(is.null(names(bias))) names(bias) <- names(x)
  stopifnot(all(names(x) %in% names(bias)))
  if(is.null(seqs)) seqs <- names(x)
  stopifnot(all(seqs %in% names(x)))
  windowSize <- as.integer(windowSize)
  
  # 1. identify candidate windows (position-wise vectors are handled one 
  # sequence at a time to limit memory usage)
  cand <- lapply(setNames(seqs,seqs), function(n){
    y <- as.integer(x[[n]])
    if(length(y) != length(bias[[n]])) 
      stop("`x` and `bias` do not have the same length for sequence ", n)
    nw <- length(y) %/% windowSize
    if(nw < 1L) return(NULL)
    y <- y[seq_len(nw*windowSize)]
    tot <- colSums(matrix(y, nrow=windowSize))
    i <- which(tot >= minCount)
    if(length(i)==0) return(NULL)
    b <- .biasToLog2(bias[[n]][seq_len(nw*windowSize)])
    mb <- colMeans(matrix(b, nrow=windowSize)[,i,drop=FALSE])
    data.frame(seq=n, win=i, meanBias=mb, stringsAsFactors=FALSE)
  })
  cand <- do.call(rbind, cand)
  if(is.null(cand) || nrow(cand) < 10) 
    stop("Not enough windows with at least ", minCount, " insertions; consider ",
         "increasing `windowSize` or lowering `minCount`.")
  
  # 2. sample windows evenly across the range of mean bias
  bin <- cut(cand$meanBias, unique(quantile(cand$meanBias, 
                                            seq(0,1,length.out=nBins+1))),
             include.lowest=TRUE, labels=FALSE)
  perBin <- ceiling(nWindows/max(bin, na.rm=TRUE))
  sel <- unlist(lapply(split(seq_len(nrow(cand)), bin), function(i){
    if(length(i) > perBin) sample(i, perBin) else i
  }), use.names=FALSE)
  cand <- cand[sel,,drop=FALSE]
  
  # 3. extract the data for the selected windows
  Y <- B <- matrix(0, nrow=windowSize, ncol=nrow(cand))
  for(n in unique(cand$seq)){
    j <- which(cand$seq==n)
    y <- as.integer(x[[n]]); b <- .biasToLog2(bias[[n]])
    idx <- outer(seq_len(windowSize), (cand$win[j]-1L)*windowSize, "+")
    Y[,j] <- y[idx]; B[,j] <- b[idx]
  }
  if(anyNA(Y) || anyNA(B)) stop("Missing values in `x` or `bias`.")
  tot <- colSums(Y)
  
  # 4. maximum conditional likelihood: the score is monotonic in k
  lw <- function(k){ # window-wise weighted mean and variance of the log2-bias
    W <- 2^(k*B); W <- t(t(W)/colSums(W))
    m <- colSums(W*B)
    list(m=m, v=colSums(W*B^2)-m^2, W=W)
  }
  sumyb <- sum(Y*B)
  score <- function(k) sumyb - sum(tot*lw(k)$m)
  if(sign(score(interval[1])) == sign(score(interval[2])))
    stop("The scale could not be estimated within the interval ", 
         paste(interval, collapse="-"))
  k <- uniroot(score, interval, tol=1e-6)$root
  info <- log(2)^2 * sum(tot * lw(k)$v)
  
  # 5. diagnostics
  b2 <- cut(as.numeric(B), unique(quantile(B, seq(0,1,length.out=nBins+1))),
            include.lowest=TRUE)
  ex <- function(k) as.numeric(t(t(lw(k)$W)*tot))
  d <- data.frame(meanBias=as.numeric(tapply(as.numeric(B), b2, mean)),
                  observed=as.numeric(tapply(as.numeric(Y), b2, sum)),
                  expected1=as.numeric(tapply(ex(1), b2, sum)),
                  expectedK=as.numeric(tapply(ex(k), b2, sum)))
  list(scale=k, se=1/sqrt(info), nWindows=ncol(Y), nInsertions=sum(tot),
       diagnostic=d)
}

# converts bias to log2(odds ratio); integers are taken to be 100*log2(OR)
.biasToLog2 <- function(b){
  if(is.integer(b)) return(as.numeric(b)/100)
  log2(as.numeric(b))
}

# reads the (shifted) insertion counts of a BAM file and the matching bias, for
# a selection of sequences
.bamInsertionsAndBias <- function(bam, genome, seqs, minCount, minMapq, 
                                  nthreads, kmerBias=NULL){
  if(is.character(bam) && grepl("\\.bam$", bam, ignore.case=TRUE)){
    if(!file.exists(paste0(bam, ".bai")) && 
       !file.exists(sub("\\.bam$", ".bai", bam, ignore.case=TRUE)))
      stop("The BAM file needs to be indexed.")
    st <- Rsamtools::idxstatsBam(bam)
    st <- st[st$mapped > 0,,drop=FALSE]
  }else if(is(bam, "GRanges")){
    tb <- table(as.character(seqnames(bam)))
    st <- data.frame(seqnames=names(tb), mapped=as.integer(tb))
  }else{ # fragment file: no read counts available
    nm <- Rsamtools::seqnamesTabix(bam)
    st <- data.frame(seqnames=nm, mapped=NA_integer_)
  }
  st$seqnames <- as.character(st$seqnames)
  gl <- seqlengths(genome)
  # map the BAM's sequence names to the genome's (e.g. "9" vs "chr9")
  gmap <- .matchSeqNames(st$seqnames, names(gl))
  if(is.null(seqs)){
    mito <- intersect(c("chrM","MT","M","chrMT"), st$seqnames)
    nm1 <- st$mapped[match(mito[1], st$seqnames)]
    if(length(mito)>0 && (is.na(nm1) || nm1 >= 100*minCount)){
      seqs <- mito[1]
    }else{
      st <- st[!is.na(gmap[st$seqnames]),,drop=FALSE]
      if(nrow(st)==0) stop("No sequence in common between the BAM and `genome`.")
      if(all(is.na(st$mapped))) st$mapped <- as.numeric(gl[gmap[st$seqnames]])
      seqs <- as.character(st$seqnames[which.max(st$mapped)])
    }
    message("Using sequence(s): ", paste(seqs, collapse=", "))
  }
  if(!all(seqs %in% st$seqnames))
    stop("Sequences without reads in the BAM file: ", 
         paste(setdiff(seqs, st$seqnames), collapse=", "))
  if(anyNA(gmap[seqs]))
    stop("Sequences not found in `genome`: ", 
         paste(seqs[is.na(gmap[seqs])], collapse=", "))
  x <- list(); bias <- list()
  for(n in seqs){
    gn <- gmap[[n]]
    L <- as.integer(gl[[gn]])
    cuts <- .getCuts(bam, GRanges(n, IRanges(1L, L)), minMapq=minMapq)
    cuts <- cuts[as.character(seqnames(cuts))==n]
    isp <- as.character(strand(cuts))=="+"
    p <- start(cuts) + ifelse(isp, 4L, -4L)
    p <- p[p >= 1L & p <= L]
    x[[n]] <- tabulate(p, nbins=L)
    rm(cuts, isp, p)
    gr <- GRanges(gn, IRanges(1L, L))
    sq <- if(is(genome, "FaFile")) Rsamtools::scanFa(genome, gr) else 
      Biostrings::getSeq(genome, gr)
    names(sq) <- n
    bias[[n]] <- getSequenceTn5Bias(sq, bias=kmerBias, nthreads=nthreads,
                                    asInteger=TRUE, verbose=FALSE)[[1]]
  }
  list(x=x, bias=bias)
}

# named vector mapping `x` to the matching name in `ref` (NA if none)
.matchSeqNames <- function(x, ref){
  cand <- function(n) c(n, paste0("chr", n), sub("^chr", "", n),
                        if(n %in% c("MT","M","chrMT")) "chrM",
                        if(n=="chrM") c("MT","M"))
  setNames(vapply(x, function(n){
    m <- intersect(cand(n), ref)
    if(length(m)==0) NA_character_ else m[1]
  }, character(1)), x)
}
