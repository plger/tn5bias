#' insertionProfiles
#'
#' Computes aggregate Tn5 cut (i.e. read 5' end) profiles around a set of 
#' regions (e.g. motif matches), both as observed and after correction for the
#' Tn5 sequence bias, separately for each strand and for each of one or more 
#' samples. The bias is computed only once per region, and shared across 
#' samples.
#' 
#' For each position of each region, the cut counts of each strand are divided
#' by the bias of that strand (see \code{\link{getSequenceTn5Bias}} with 
#' `shift=FALSE`), raised to the sample's scale (see 
#' \code{\link{estimateBiasScale}}), and then summed across regions. Note that 
#' the bias, and therefore the correction, is applied on each cut separately 
#' and is thus unaffected by the regions' orientation. Strands are reported 
#' relative to the orientation of the regions: for regions on the negative 
#' strand, the profile is reversed and the strands are swapped, so that 
#' strand "+" always means "same strand as the region".
#'
#' @param regions A \code{\link[GenomicRanges]{GRanges}} of regions (e.g. motif 
#'   matches). Each is resized to `2*flank+1` around its center. Regions that 
#'   extend beyond the sequences' boundaries are dropped.
#' @param samples A named character vector of paths to the samples' files, 
#'   which can be either indexed BAM files (the 5' ends of all reads are used, 
#'   so it works for single- as well as paired-end data) or tabix-indexed 
#'   fragment files (ending with .tsv.gz or .bed.gz). Alternatively, a named 
#'   list of per-sample \code{GRanges} of width-1 cut sites, with the strand 
#'   indicating which strand the cut is on.
#' @param genome A \link[Rsamtools]{FaFile} or \code{BSgenome} from which to 
#'   fetch the sequences. Sequence names must match those of `regions`.
#' @param flank The number of nucleotides on each side of the region centers 
#'   (default 100).
#' @param scale The bias scale, either a single number or a (named) vector 
#'   with one value per sample. See \code{\link{estimateBiasScale}}.
#' @param bias An optional 10-mer bias table (see \code{\link{getTn5Bias}}).
#' @param minMapq For BAM files, the minimum mapping quality (default 10).
#' @param barcodes For fragment files, an optional vector of cell barcodes to 
#'   include. 
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list with the following elements:
#' \describe{
#'  \item{`raw` and `corrected`}{ Arrays of dimensions positions x strand x 
#'   sample with the aggregate observed and bias-corrected cut counts.}
#'  \item{`bias`}{ A matrix (positions x strand) with the mean bias (as odds
#'   ratio, with scale 1) across regions.}
#'  \item{`nRegions`}{ The number of regions used (after dropping those 
#'   out-of-bounds).}
#'  \item{`scale`}{ The scale used for each sample.}
#'  \item{`positions`}{ The position relative to the region center.}
#' }
#' @export
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlengths seqnames
#' @importFrom BiocGenerics strand start end width
#' @importFrom Rsamtools scanFa scanBam ScanBamParam scanBamFlag scanTabix
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom utils read.delim
insertionProfiles <- function(regions, samples, genome, flank=100L, scale=1, 
                              bias=NULL, minMapq=10L, barcodes=NULL, 
                              verbose=TRUE){
  stopifnot(is(regions, "GRanges"), length(regions) > 0, is.numeric(scale))
  flank <- as.integer(flank)
  if(is.null(names(samples))) names(samples) <- 
    if(is.character(samples)) basename(samples) else paste0("sample", seq_along(samples))
  if(length(scale)==1) scale <- setNames(rep(scale, length(samples)), names(samples))
  if(is.null(names(scale))) names(scale) <- names(samples)
  stopifnot(setequal(names(scale), names(samples)))
  scale <- scale[names(samples)]
  pad <- 20L # flank needed to compute the bias at the edges
  W <- 2L*flank + 1L
  
  # center and resize regions, drop those out of bounds
  isneg <- as.character(strand(regions))=="-"
  ctr <- as.integer(floor((start(regions) + end(regions))/2))
  sl <- if(is(genome, "FaFile")) seqlengths(genome) else seqlengths(genome)
  sn <- as.character(seqnames(regions))
  ok <- sn %in% names(sl) & ctr-flank-pad >= 1L
  ok[ok] <- ctr[ok]+flank+pad <= sl[sn[ok]]
  if(!any(ok)) stop("No valid region.")
  if(verbose && !all(ok)) message(sum(!ok), " regions were dropped (out of bounds).")
  sn <- sn[ok]; ctr <- ctr[ok]; isneg <- isneg[ok]
  nr <- length(ctr)
  gr <- GRanges(sn, IRanges(ctr-flank, ctr+flank))
  grx <- GRanges(sn, IRanges(ctr-flank-pad, ctr+flank+pad))
  
  # bias, computed once for all samples
  if(verbose) message("Computing bias for ", nr, " regions...")
  sq <- if(is(genome, "FaFile")) scanFa(genome, grx) else 
    Biostrings::getSeq(genome, grx)
  bm <- function(complement){
    b <- getSequenceTn5Bias(sq, bias=bias, shift=FALSE, complement=complement,
                            verbose=FALSE)
    matrix(unlist(b, use.names=FALSE), nrow=W+2L*pad)[pad+seq_len(W),,drop=FALSE]
  }
  bias_g <- list(plus=bm(FALSE), minus=bm(TRUE)) # in genomic orientation
  rm(sq)
  
  # orientation: for negative-strand regions, reverse and swap strands
  orient <- function(plus, minus){
    rv <- function(m) if(any(isneg)) {m[,isneg] <- m[rev(seq_len(W)),isneg]; m} else m
    plus <- rv(plus); minus <- rv(minus)
    if(any(isneg)){
      tmp <- plus[,isneg]; plus[,isneg] <- minus[,isneg]; minus[,isneg] <- tmp
    }
    list(plus=plus, minus=minus)
  }
  b <- orient(bias_g$plus, bias_g$minus)
  rm(bias_g)
  
  strands <- c("+","-")
  pos <- seq(-flank, flank)
  dn <- list(position=pos, strand=strands, sample=names(samples))
  raw <- corr <- array(0, dim=c(W, 2L, length(samples)), dimnames=dn)
  
  for(s in names(samples)){
    if(verbose) message("Counting cuts in ", s, "...")
    cuts <- .getCuts(samples[[s]], grx, minMapq=minMapq, barcodes=barcodes)
    cn <- .countCuts(cuts, gr)  # list(plus=, minus=) in genomic orientation
    rm(cuts)
    cn <- orient(cn$plus, cn$minus)
    for(i in 1:2){
      st <- c("plus","minus")[i]
      raw[,i,s] <- rowSums(cn[[st]])
      corr[,i,s] <- rowSums(cn[[st]] / b[[st]]^scale[[s]])
    }
  }
  list(raw=raw, corrected=corr, 
       bias=cbind("+"=rowMeans(b$plus), "-"=rowMeans(b$minus)),
       nRegions=nr, scale=scale, positions=pos)
}

# returns a GRanges of width-1 cut sites (strand = strand of the cut)
.getCuts <- function(x, regions, minMapq=10L, barcodes=NULL){
  if(is(x, "GRanges")) return(x)
  stopifnot(is.character(x), length(x)==1)
  rr <- reduce(regions)
  if(grepl("\\.bam$", x, ignore.case=TRUE)){
    ga <- readGAlignments(x, param=ScanBamParam(which=rr, 
                   flag=scanBamFlag(isUnmappedQuery=FALSE, 
                                    isSecondaryAlignment=FALSE,
                                    isDuplicate=FALSE),
                   mapqFilter=minMapq))
    isp <- as.character(strand(ga))=="+"
    sn <- as.character(seqnames(ga))
    p <- ifelse(isp, start(ga), end(ga))
    return(GRanges(sn, IRanges(p, width=1L), 
                   strand=ifelse(isp, "+", "-")))
  }
  # fragment file
  lines <- unlist(scanTabix(x, param=rr), use.names=FALSE)
  if(length(lines)==0) return(GRanges())
  d <- read.delim(text=lines, header=FALSE, stringsAsFactors=FALSE,
                  colClasses=c("character","integer","integer",
                               rep("NULL", 0)), 
                  nrows=length(lines))
  if(!is.null(barcodes)){
    bc <- vapply(strsplit(lines, "\t", fixed=TRUE), `[`, character(1), 4L)
    d <- d[bc %in% barcodes,,drop=FALSE]
  }
  # fragment files are 0-based, end-exclusive
  c(GRanges(d[[1]], IRanges(d[[2]]+1L, width=1L), strand="+"),
    GRanges(d[[1]], IRanges(d[[3]], width=1L), strand="-"))
}

# counts the cuts at each position of (same-width) regions, per strand
.countCuts <- function(cuts, regions){
  W <- width(regions)[1]
  res <- lapply(c("+","-"), function(s){
    cs <- cuts[as.character(strand(cuts))==s]
    if(length(cs)==0) return(matrix(0L, nrow=W, ncol=length(regions)))
    ov <- findOverlaps(cs, regions)
    # position within region, region index
    ri <- subjectHits(ov)
    pi <- start(cs)[queryHits(ov)] - start(regions)[ri] + 1L

    k <- tabulate((ri-1L)*W + pi, nbins=W*length(regions))
    matrix(k, nrow=W)
  })
  list(plus=res[[1]], minus=res[[2]])
}
