## We need to check that insertions and deletions work separately and together
## If possible, check the orders of the returned GRanges
## No sequence analysis is really needed so just create some dummy ranges

var <- c(
    GRanges("seq1:6-10"), # Deletion
    GRanges("seq1:12") #Insertion
)
var$REF <- c("AGCTA", "C")
var$ALT <- c("A", "CCC")
var$change <- nchar(var$ALT) - nchar(var$REF)
var$shift <- cumsum(var$change)
gr <- c(
    GRanges("seq1:1-10:+"),
    GRanges("seq1:11:+"),
    GRanges("seq1:21-25:-"),
    GRanges("seq1:1-20:*")
)
gr$feature <- paste0("feature", seq_along(gr))

test_that("Errors are as expected", {
    expect_error(shiftByVar(gr, var), "The seqinfo.+")
})

## Now add the sequenq lengths
seqlengths(gr) <- c(seq1 = 50)

test_that("Shifts are numerically correct", {
    new_gr <- shiftByVar(gr, var)
    ## Step through each range checking for correctness
    expect_equal(width(new_gr[1]), width(gr[1]) + var$change[1])
    expect_equal(start(new_gr[2]), start(gr[2]) + var$change[1])
    expect_equal(width(new_gr)[2:3], width(gr)[2:3])
    expect_equal(width(new_gr[4]), width(gr)[4] + sum(var$change))
    expect_equal(seqlengths(new_gr), seqlengths(gr) + sum(var$change))
})
