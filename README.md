# tn5bias : Computing Tn5 bias across sequences for ATAC-seq correction

** Very much under development, do not use! **

ATAC-seq (Assay for Transposase-Accessible Chromatin) maps accessible regions of
the genome through the use of a hyperactive transposase, Tn5, which inserts 
itself in (accessible regions of) the genome, and in doing so inserts sequencing
primers on each side of the insertion sites. The resulting fragments can be 
sequenced, providing an idea of the degree of accessibility of different regions.

However, the Tn5 itself has some (weak) sequence bias beyond accessibility. The 
`tn5bias` R package uses kmer bias estimated by SELMA (Hu et al., 2022) to 
compute a bias for each position in a given DNA sequence, taking both strands 
into account. Such a bias profile can then be used to correct ATAC-seq profiles.
