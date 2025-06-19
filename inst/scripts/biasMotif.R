data(tn5bias, package="tn5bias")

# build kmers:
nts <- c("A","C","G","T")
kmers <- expand.grid(nts,nts,nts,nts,nts,nts,nts,nts,nts,nts)
bias <- tn5bias/100
tn5biasMotif <- sapply(seq_len(10), \(j){
  rowsum(bias, kmers[,j])/262144
})
row.names(tn5biasMotif) <- nts
#save(tn5biasMotif, file="../../data/tn5biasMotif.RData")
