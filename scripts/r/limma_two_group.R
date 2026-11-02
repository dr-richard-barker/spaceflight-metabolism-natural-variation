# limma moderated t for one two-group comparison on log2 expression (RMA microarray values or protein
# abundances), for prereg/revision1 R4 and R5. Called by scripts/26 and scripts/27:
#   Rscript limma_two_group.R expr.tsv design.tsv out.tsv
# expr.tsv: features x samples; design.tsv: columns sample, group (treated/control).
# Writes, per feature, logFC, moderated t and p for treated - control.
suppressPackageStartupMessages(library(limma))
a <- commandArgs(trailingOnly = TRUE)
x <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
des <- read.delim(a[2], stringsAsFactors = FALSE)
stopifnot(identical(des$sample, colnames(x)), setequal(unique(des$group), c("treated", "control")))
group <- factor(des$group, levels = c("control", "treated"))
design <- model.matrix(~ group)
fit <- eBayes(lmFit(x, design))
tt <- topTable(fit, coef = "grouptreated", number = Inf, sort.by = "none")
out <- data.frame(feature = rownames(tt), logFC = tt$logFC, t = tt$t, P = tt$P.Value, check.names = FALSE)
write.table(out, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  limma %s, %d features, %d v %d samples, df.prior %.2f\n", packageVersion("limma"), nrow(out),
            sum(group == "treated"), sum(group == "control"), fit$df.prior))
