#' @title Tn5 bias
#' @description
#' Raw Tn5 bias per 10-mers obtained using SELMA (Hu et al. 2022),
#' downloaded from the SELMA repository. For low memory footprint
#' this has been saved as integers without the kmers, using
#' system.file('inst/scripts/', 'bias_data_gen.R', package='tn5bias')
#' Use \link{getTn5Bias} to obtain an interpretable kmer-bias table.
"tn5bias"

#' @title Tn5 bias motif
#' @description
#' A motif summary of the Tn5 bias, summarized from the kmer-level bias
#' obtained using SELMA (Hu et al. 2022). See \link{getTn5Bias} for more
#' information. This has been generated using
#' system.file('inst/scripts/', 'biasMotif.R', package='tn5bias')
"tn5biasMotif"

