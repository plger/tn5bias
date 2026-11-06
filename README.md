# tn5bias : Computing Tn5 bias across sequences for ATAC-seq correction

ATAC-seq (Assay for Transposase-Accessible Chromatin) maps accessible regions of
the genome through the use of a hyperactive transposase, Tn5, which inserts 
itself in (accessible regions of) the genome, and in doing so inserts sequencing
primers on each side of the insertion sites. The resulting fragments can be 
sequenced, providing an idea of the degree of accessibility of different regions.

However, the Tn5 itself has some sequence bias beyond accessibility. The 
`tn5bias` R package uses kmer bias estimated by SELMA (Hu et al., 2022) to 
compute a bias for each position in a given DNA sequence, taking both strands 
into account. Such bias can then be scaled to the magnitude it has in given 
samples. The assumption (supported by the literature), is that the bias is a 
property of the enzyme which only changes in magnitude across samples or 
experiments.

## Getting started

Install using

```
BiocManager::install("plger/tn5bias")
```

Then, assuming a sequence `seq`, compute its nucleotide-level bias using:

```
bias <- getSequenceTn5Bias(seq)
```

Such bias can then be used to correct ATAC-seq profiles.

For more detail, and functions related to specific use-cases, see the vignette:

```
vignette("tn5bias")
```
