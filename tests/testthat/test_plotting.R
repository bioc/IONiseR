library(IONiseR)
context('Testing plotting functions')

require(minionSummaryData)
data(s.typhi.rep2, package = "minionSummaryData")

p1 <- plotReadCategoryCounts(s.typhi.rep2)
test_that("plot is ggplot", {
    expect_is(p1, "ggplot")
})

test_that("category counts identified", {
    expect_equal(p1$data[['count']], c(3591, 3088, 1678, 1412, 1019))
    expect_equal(p1$data[['category']], 
                 factor(c("Fast5 File Count", "Template", "Complement", "Full 2D", "Pass"),
                        levels = c("Fast5 File Count", "Template", "Complement", "Full 2D", "Pass")))
})

p2 <- plotActiveChannels(s.typhi.rep2)
test_that("plot is ggplot", {
    expect_is(p2, "ggplot")
})

p3 <- plotReadAccumulation(s.typhi.rep2)
test_that("plot is ggplot", {
    expect_is(p3, "ggplot")
})

p4 <- plotReadTypeProduction(s.typhi.rep2)
test_that("plot is ggplot", {
    expect_is(p4, "ggplot")
})

p5 <- plotKmerFrequencyCorrelation(s.typhi.rep2)
test_that("plot is ggplot", {
    expect_is(p5, "ggplot")
})