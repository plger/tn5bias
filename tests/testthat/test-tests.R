bias <- getTn5Bias()

test_that("kmer table construction works", {
  expect_false(any(is.na(bias$bias)))
})

test_that("bias estimation works", {
  seq <- "ACGTCGCGACGTTTGCATATACGGCGCA"
  seqbias <- getSequenceTn5Bias(seq, bias=bias, verbose=FALSE)
  expect_true(is(seqbias, "NumericList"))
  expect_false(any(is.na(seqbias[[1]])))
  expect_equal(length(seqbias[[1]]), nchar(seq))
})

test_that("some extra params works", {
  seq <- Biostrings::DNAString("ACGTCGCGACGTTTGCATATACGGCGCA")
  seqbias <- getSequenceTn5Bias(seq, bias=bias, shift=FALSE,
                                complement=TRUE, asInteger=TRUE)
  expect_true(is(seqbias, "IntegerList"))
  expect_false(any(is.na(seqbias[[1]])))
  expect_equal(length(seqbias[[1]]), length(seq))
})

test_that("bias matches a direct implementation of SELMA", {
  set.seed(1)
  s <- paste(sample(c("A","C","G","T"), 200, replace=TRUE), collapse="")
  lb <- setNames(log2(bias$bias), bias$kmer)
  rc <- function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))
  sb <- function(a, l=10) substr(s, a, a+l-1)
  # cuts at 1-based position p (g=p-1 0-based)
  plus <- function(p) sqrt(2^lb[sb(p-5)] * 2^lb[rc(sb(p+4))])
  minus <- function(p) sqrt(2^lb[rc(sb(p-4))] * 2^lb[sb(p-13)])
  pos <- 30:150
  f <- function(...) as.numeric(getSequenceTn5Bias(s, bias=bias, verbose=FALSE, ...)[[1]])
  expect_equal(f(shift=FALSE)[pos], unname(sapply(pos, plus)))
  expect_equal(f(shift=FALSE, complement=TRUE)[pos], unname(sapply(pos, minus)))
  expect_equal(f()[pos], unname(sapply(pos, function(p) plus(p-4))))
  # both strands are the same once shifted
  expect_equal(f()[pos], unname(sapply(pos, function(p) minus(p+4))))
})

test_that("FaFile input works", {
  set.seed(2)
  seqs <- Biostrings::DNAStringSet(c(a=paste(sample(c("A","C","G","T"), 100, TRUE), collapse=""),
                                     b=paste(sample(c("A","C","G","T"), 60, TRUE), collapse="")))
  f <- tempfile(fileext=".fa")
  Biostrings::writeXStringSet(seqs, f)
  Rsamtools::indexFa(f)
  expect_equal(getSequenceTn5Bias(Rsamtools::FaFile(f), bias=bias, verbose=FALSE),
               getSequenceTn5Bias(seqs, bias=bias, verbose=FALSE))
})

test_that("scale amplifies the bias and can be estimated", {
  set.seed(3)
  s <- paste(sample(c("A","C","G","T"), 30000, replace=TRUE), collapse="")
  b1 <- getSequenceTn5Bias(s, bias=bias, verbose=FALSE)
  b2 <- getSequenceTn5Bias(s, bias=bias, verbose=FALSE, scale=1.7)
  expect_equal(log2(as.numeric(b2[[1]])), 1.7*log2(as.numeric(b1[[1]])))
  p <- log2(as.numeric(b2[[1]]))
  y <- rpois(length(p), 0.5 * 2^p / mean(2^p))
  est <- estimateBiasScale(list(seq1=y), b1, minCount=10)
  expect_lt(abs(est$scale - 1.7), 0.15)
  # integer-encoded bias gives the same result
  b1i <- getSequenceTn5Bias(s, bias=bias, verbose=FALSE, asInteger=TRUE)
  expect_lt(abs(estimateBiasScale(list(seq1=y), b1i, minCount=10)$scale - est$scale), 0.05)
  expect_error(estimateBiasScale(list(seq1=y[1:100]), b1i[1], minCount=10))
})

test_that("insertionProfiles works across strands, samples and input types", {
  library(GenomicRanges)
  set.seed(4)
  chr <- paste(sample(c("A","C","G","T"), 5000, replace=TRUE), collapse="")
  f <- tempfile(fileext=".fa")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(chrT=chr)), f)
  Rsamtools::indexFa(f); fa <- Rsamtools::FaFile(f)
  mkcuts <- function(n){
    GRanges("chrT", IRanges(sample(1000:4000, n, replace=TRUE), width=1),
            strand=sample(c("+","-"), n, replace=TRUE))
  }
  s1 <- mkcuts(20000); s2 <- mkcuts(5000)
  regs <- GRanges("chrT", IRanges(c(1500, 2500, 3500), width=1), strand="+")
  r <- insertionProfiles(regs, list(a=s1, b=s2), fa, flank=50, 
                         scale=c(a=1, b=1.5), bias=bias, verbose=FALSE)
  expect_equal(dim(r$raw), c(101, 2, 2))
  expect_equal(dimnames(r$raw)$sample, c("a","b"))
  expect_equal(r$nRegions, 3)
  # raw counts are consistent with the cuts
  tot <- sum(r$raw[,,"a"])
  expect_equal(tot, sum(overlapsAny(s1, resize(regs, 101, fix="center"), ignore.strand=TRUE)))
  # corrected = raw/bias for a single region
  r1 <- insertionProfiles(regs[1], list(a=s1), fa, flank=50, bias=bias, verbose=FALSE)
  bp <- as.numeric(getSequenceTn5Bias(chr, bias=bias, shift=FALSE, verbose=FALSE)[[1]])
  rawp <- r1$raw[,"+","a"]
  expect_equal(unname(r1$corrected[,"+","a"]), unname(rawp/bp[1450:1550]))
  # every region has its own bias (sequence names are duplicated)
  r3 <- insertionProfiles(regs[3], list(a=s1), fa, flank=50, bias=bias, verbose=FALSE)
  expect_equal(unname(r3$corrected[,"+","a"]), 
               unname(r3$raw[,"+","a"]/bp[3450:3550]))
  expect_equal(unname(r$bias[,"+"]), 
               unname(rowMeans(sapply(c(1500,2500,3500), function(p) bp[(p-50):(p+50)]))))
  # negative-strand region: reversed and strand-swapped
  rneg <- insertionProfiles(regs[1], list(a=s1), fa, flank=50, bias=bias, verbose=FALSE)
  regn <- regs[1]; strand(regn) <- "-"
  rn <- insertionProfiles(regn, list(a=s1), fa, flank=50, bias=bias, verbose=FALSE)
  expect_equal(unname(rn$raw[,"+","a"]), rev(unname(rneg$raw[,"-","a"])))
  expect_equal(unname(rn$corrected[,"-","a"]), rev(unname(rneg$corrected[,"+","a"])))
  # fragment file input is equivalent
  fr <- tempfile(fileext=".bed")
  set.seed(5)
  st <- sample(1000:3900, 3000, replace=TRUE); en <- st + sample(40:200, 3000, replace=TRUE)
  d <- data.frame("chrT", st-1L, en, "bc", 1L)
  d <- d[order(d[,2]),]
  write.table(d, fr, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  gz <- Rsamtools::bgzip(fr, overwrite=TRUE); Rsamtools::indexTabix(gz, format="bed")
  cuts <- c(GRanges("chrT", IRanges(st, width=1), strand="+"),
            GRanges("chrT", IRanges(en, width=1), strand="-"))
  ra <- insertionProfiles(regs, list(a=gz), fa, flank=50, bias=bias, verbose=FALSE)
  rb <- insertionProfiles(regs, list(a=cuts), fa, flank=50, bias=bias, verbose=FALSE)
  expect_equal(ra$raw, rb$raw)
  # out-of-bounds regions are dropped
  expect_equal(insertionProfiles(c(regs, GRanges("chrT", IRanges(10, width=1))), 
                                 list(a=cuts), fa, flank=50, bias=bias,
                                 verbose=FALSE)$nRegions, 3)
})

test_that("getRegionTn5Bias works", {
  library(GenomicRanges)
  set.seed(6)
  chr <- paste(sample(c("A","C","G","T"), 3000, replace=TRUE), collapse="")
  f <- tempfile(fileext=".fa")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(chrT=chr)), f)
  Rsamtools::indexFa(f); fa <- Rsamtools::FaFile(f)
  regs <- GRanges(c("chrT","chrT","chrT","chrX"), 
                  IRanges(c(500, 1000, 2950, 10), width=c(100, 300, 200, 50)))
  b <- as.numeric(getSequenceTn5Bias(chr, bias=bias, verbose=FALSE)[[1]])
  r <- getRegionTn5Bias(regs, fa, bias=bias, verbose=FALSE)
  expect_equal(r[1], mean(b[500:599]))
  expect_equal(r[2], mean(b[1000:1299]))
  # edge: positions without a full flank are ignored (no padding bias of 1)
  r3 <- getRegionTn5Bias(regs[3], fa, bias=bias, maxNFraction=1, verbose=FALSE)
  expect_equal(r3, mean(b[2950:2991]))
  expect_true(is.na(r[3]))  # 75% of the region is outside the sequence
  expect_true(is.na(r[4]))
  expect_equal(getRegionTn5Bias(regs[1:2], fa, bias=bias, method="geometric", 
                                log2=TRUE, verbose=FALSE),
               c(mean(log2(b[500:599])), mean(log2(b[1000:1299]))))
  expect_equal(getRegionTn5Bias(regs[1], fa, bias=bias, scale=2, 
                                method="geometric", log2=TRUE, verbose=FALSE),
               2*mean(log2(b[500:599])))
})

test_that("getRegionTn5Bias supports several scores", {
  set.seed(7)
  chr <- paste(sample(c("A","C","G","T"), 2000, replace=TRUE), collapse="")
  f <- tempfile(fileext=".fa")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(chrT=chr)), f)
  Rsamtools::indexFa(f); fa <- Rsamtools::FaFile(f)
  regs <- GenomicRanges::GRanges("chrT", IRanges::IRanges(c(300, 900), width=300))
  b <- as.numeric(getSequenceTn5Bias(chr, bias=bias, verbose=FALSE)[[1]])
  r <- getRegionTn5Bias(regs, fa, bias=bias, method=c("mean","max","quantile"), 
                        q=c(0.9,0.98), verbose=FALSE)
  expect_equal(colnames(r), c("mean","max","quantile_0.9","quantile_0.98"))
  expect_equal(unname(r[1,"max"]), max(b[300:599]))
  expect_equal(unname(r[2,"quantile_0.98"]), unname(quantile(b[900:1199], 0.98)))
  expect_equal(unname(r[,"mean"]), 
               getRegionTn5Bias(regs, fa, bias=bias, verbose=FALSE))
})

test_that("unknown bases are handled", {
  set.seed(8)
  chr <- paste(sample(c("A","C","G","T"), 2000, replace=TRUE), collapse="")
  substr(chr, 1000, 1049) <- strrep("N", 50)
  f <- tempfile(fileext=".fa")
  Biostrings::writeXStringSet(Biostrings::DNAStringSet(c(chrT=chr)), f)
  Rsamtools::indexFa(f); fa <- Rsamtools::FaFile(f)
  regs <- GenomicRanges::GRanges("chrT", IRanges::IRanges(c(900, 1000, 1500, 1020), 
                                                          width=c(300, 100, 200, 20)))
  b <- as.numeric(getSequenceTn5Bias(chr, bias=bias, verbose=FALSE)[[1]])
  r <- getRegionTn5Bias(regs, fa, bias=bias, verbose=FALSE)
  expect_true(all(is.na(r[c(2,4)])))            # (almost) all N
  expect_false(is.na(r[3]))
  # 900-1199: N at 1000-1049 (n=50, 17%): not scored with the default, 
  expect_true(is.na(r[1]))
  # ... but scored if allowed, excluding positions within 9bp of an N
  r2 <- getRegionTn5Bias(regs[1], fa, bias=bias, maxNFraction=0.5, verbose=FALSE)
  valid <- setdiff(900:1199, 991:1058)
  expect_equal(r2, mean(b[valid]))
  expect_equal(r[3], mean(b[1500:1699]))
})

test_that("scaleBias is equivalent to scaling at computation", {
  s <- "ACGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGATCGATTTGACCGATAGCG"
  b <- getSequenceTn5Bias(s, bias=bias, verbose=FALSE)
  expect_equal(scaleBias(b, 1.5), 
               getSequenceTn5Bias(s, bias=bias, scale=1.5, verbose=FALSE))
  bi <- getSequenceTn5Bias(s, bias=bias, verbose=FALSE, asInteger=TRUE)
  expect_true(is(scaleBias(bi, 1.5), "IntegerList"))
  expect_equal(as.numeric(scaleBias(bi, 1.5)[[1]]), 
               as.numeric(getSequenceTn5Bias(s, bias=bias, scale=1.5, verbose=FALSE,
                                             asInteger=TRUE)[[1]]), tolerance=1/100*2)
  expect_equal(scaleBias(c(1,2,4), 2), c(1,4,16))
  expect_equal(scaleBias(c(100L,-50L), 1.5), c(150L,-75L))
})

test_that("duplicated sequence names are handled", {
  x <- Biostrings::DNAStringSet(c(a="ACGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGAT", 
                                  a="TTGACCGATAGCGCGCTATATCGGCGATATCGCGATA"))
  b <- getSequenceTn5Bias(x, bias=bias, verbose=FALSE)
  b2 <- getSequenceTn5Bias(x[2], bias=bias, verbose=FALSE)
  expect_equal(as.numeric(b[[2]]), as.numeric(b2[[1]]))
})
