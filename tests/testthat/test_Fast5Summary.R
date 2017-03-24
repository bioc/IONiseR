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

test_that("subsetting works", {
    expect_equal(length(s.typhi.rep2[1:10]), 10)
})

fast5files <- system.file('extdata', c('example.fast5', 'example_v2.fast5'), package = "IONiseR")

test_that("Error thrown when incompatible files combined", {
    expect_error(readFast5Summary(fast5files), "Inconsistent analysis workflows")
})

test_that("Error thrown when file can't be found", {
    expect_error(readFast5Summary('/not/path/to/file.fast5'), "None of the provided files can be accessed")
})

f1 <- readFast5Summary(fast5files[2])

test_that("Can fast5 file be read?", {
  expect_is(f1, 'Fast5Summary')
})

test_that("accessors", {
    expect_is(readInfo(f1), 'data.frame')
    expect_is(eventData(f1), 'data.frame')
    expect_true(all(dim(eventData(f1)) == c(1,4)))
    expect_is(baseCalled(f1), 'data.frame')
    expect_is(fastq(f1), 'ShortReadQ')
})

test_that("Show method prints summary", {
    expect_output(show(f1), regexp = "^Object of class: Fast5Summary")
})

test_that("Reading FASTQ", {
    expect_equal(length(fastq2D(f1)), 1)
    expect_equal(length(fastqTemplate(f1)), 1)
    expect_equal(length(fastqComplement(f1)), 1)
})

test_that("Catch broken file", {
    fast5_nw <- system.file('extdata', 'example_not-working.fast5', package = "IONiseR")
    expect_error(readFast5Summary(fast5_nw), "No basecalls for template strand found")
})