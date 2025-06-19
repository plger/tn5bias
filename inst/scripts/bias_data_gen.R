tmp <- tempfile(fileext = ".txt.gz")
download.file("https://github.com/zang-lab/SELMA/raw/refs/heads/master/lib/refdata/ATAC_SELMAbias_10mer.txt.gz", tmp)
b <- read.delim(tmp, header=FALSE, row.names = 1)

# build kmers:
nts <- c("A","C","G","T")
kmers <- expand.grid(nts,nts,nts,nts,nts,nts,nts,nts,nts,nts)
kmersT <- do.call(paste0, as.list(as.data.frame(kmers)))

# order bias data in the built order:
b <- b[kmersT,]
stopifnot(all(!is.na(b)))
b <- as.integer(round(100*b))

kmers2 <- unlist(lapply(1:5, \(i){
  k2 <- do.call(expand.grid, lapply((i+1):10, \(x) nts))
  Nstring <- paste0(rep("N",i), collapse="")
  k2 <- paste0(Nstring, do.call(paste0, as.list(as.data.frame(k2))))
  ag <- aggregate(b, by=kmers[,-seq_len(i)], FUN=mean)
  row.names(ag) <- paste0(Nstring, do.call(paste0, as.list(ag[,-ncol(ag)])))
  as.integer(round(ag[k2,ncol(ag)]))
}))
kmers3 <- unlist(lapply(6:10, \(i){
  k2 <- do.call(expand.grid, lapply(1:(i-1), \(x) nts))
  Nstring <- paste0(rep("N",i-5L), collapse="")
  k2 <- paste0(do.call(paste0, as.list(as.data.frame(k2))), Nstring)
  ag <- aggregate(b, by=kmers[,1:(i-1)], FUN=mean)
  row.names(ag) <- paste0(do.call(paste0, as.list(ag[,1:(i-1)])), Nstring)
  as.integer(round(ag[k2,ncol(ag)]))
}))

tn5bias <- c(b, kmers2, kmers3) 

#save(tn5bias, file="../../data/tn5bias.RData")
