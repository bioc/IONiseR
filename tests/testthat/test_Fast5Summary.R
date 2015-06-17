library(IONiseR)
context('Testing methods for the Fast5Summary class')

require(minionSummaryData)
data(s.typhi.rep2, package = "minionSummaryData")

test_that("length is the number of files", {
  expect_equal(length(s.typhi.rep2), 3591)
})

test_that("dim is not defined", {
  expect_equal(dim(s.typhi.rep2), NULL)
})

fast5file <- system.file('extdata', 'example.fast5', package = "IONiseR")
test_that("Can fast5 file be read?", {
  expect_is(readFast5Summary(fast5file), 'Fast5Summary')
})