
## The first read here is missing a quality value, so is invalid
testString <- c("read1\nAAAAAA\n+\n!!!!!", "read2\nCCCCCC\n+\n!!!!!!", NULL)

fastqResult <- .processFastqVec(testString)
expect_equal(length(fastqResult$invalid), 1)

