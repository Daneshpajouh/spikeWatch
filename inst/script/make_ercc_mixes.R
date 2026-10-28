# Builds inst/extdata/ercc_mixes.tsv, the ERCC ExFold Mix 1 and Mix 2
# design used by erccDesign().
#
# Source: the data object ERCCMix1and2 in the Bioconductor package
# erccdashboard (licence GPL (>= 2)), read from version 1.46.0. It lists
# 96 ERCC sequences with their Mix 1 and Mix 2 concentrations in
# attomoles per microlitre and the Mix 1 to Mix 2 ratio subpool. Four
# sequences have no subpool and a concentration of 0; they are not in the
# mixes and are dropped, leaving 92.
#
# Checks: 92 rows, 23 per subpool, and every Mix 1 / Mix 2 ratio equal to
# its subpool label within a relative 1e-5. The concentrations are rounded
# in the source, so the largest relative difference is about 5e-7.
library(erccdashboard)
data(ERCC)
t <- ERCCMix1and2
t <- t[!is.na(t$Subpool) & t$Mix1Conc.Attomoles_ul > 0, ]
out <- data.frame(feature = as.character(t$ERCC.AMB.Expected),
                  subgroup = as.character(t$Subpool),
                  Mix1 = t$Mix1Conc.Attomoles_ul,
                  Mix2 = t$Mix2Conc.Attomoles_ul)
stopifnot(nrow(out) == 92, all(table(out$subgroup) == 23))
label <- c(`4:1` = 4, `1:1` = 1, `1:1.5` = 1 / 1.5, `1:2` = 0.5)
stopifnot(isTRUE(all.equal(out$Mix1 / out$Mix2,
                           unname(label[out$subgroup]),
                           tolerance = 1e-5)))
write.table(out, "inst/extdata/ercc_mixes.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
