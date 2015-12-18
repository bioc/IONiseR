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

fast5files <- system.file('extdata', c('example.fast5', 'example_v2.fast5'), package = "IONiseR")
f1 <- readFast5Summary(fast5files)

test_that("Can fast5 file be read?", {
  expect_is(f1, 'Fast5Summary')
})

test_that("acessors work", {
    expect_is(readInfo(f1), 'data.frame')
    expect_is(rawData(f1), 'data.frame')
    expect_true(all(dim(rawData(f1)) == c(2,5)))
    expect_is(baseCalled(f1), 'data.frame')
    expect_is(fastq(f1), 'ShortReadQ')
})

test_that("Show method prints summary", {
    expect_output(show(f1), regexp = "^Object of class: Fast5Summary")
})