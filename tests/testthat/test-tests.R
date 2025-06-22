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

mo <- matrix( c(0,0,0,100,
                100,0,0,0,
                0,60,40,0,
                0,100,0,0,
                0,0,100,0,
                0,50,50,0), nrow=4)
row.names(mo) <- c("A","C","G","T")

  test_that("motif bias estimation works", {
  msb <- getMotifTn5Bias(mo, bias=bias, verbose=FALSE)
  expect_true(is.numeric(msb))
  expect_false(any(is.na(msb)))
  expect_equal(length(msb), ncol(mo))
})

test_that("some extra params works", {
  seq <- Biostrings::DNAString("ACGTCGCGACGTTTGCATATACGGCGCA")
  seqbias <- getSequenceTn5Bias(seq, bias=bias, shift=FALSE,
                                complement=TRUE, asInteger=TRUE)
  expect_true(is(seqbias, "IntegerList"))
  expect_false(any(is.na(seqbias[[1]])))
  expect_equal(length(seqbias[[1]]), length(seq))
})
