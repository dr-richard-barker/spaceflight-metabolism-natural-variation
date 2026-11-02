# limma-voom with edgeR TMM on a two-factor design (flight x light), for prereg/revision1 R3.
# Called by scripts/24_osd678_interaction.py:
#   Rscript limma_voom_contrasts.R counts.tsv design.tsv out.tsv
# counts.tsv: genes x samples; design.tsv: columns sample, flight (flight/ground), light (light/dark).
# Writes, per gene, logFC, moderated t and p for flight-vs-ground in light, in dark, and their
# difference (the interaction).
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
})
a <- commandArgs(trailingOnly = TRUE)
counts <- as.matrix(read.delim(a[1], row.names = 1, check.names = FALSE))
des <- read.delim(a[2], stringsAsFactors = FALSE)
stopifnot(identical(des$sample, colnames(counts)))
group <- factor(paste(des$flight, des$light, sep = "_"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
y <- normLibSizes(DGEList(counts = counts), method = "TMM")
v <- voom(y, design)
fit <- lmFit(v, design)
cm <- makeContrasts(flight_light = flight_light - ground_light,
                    flight_dark = flight_dark - ground_dark,
                    interaction = (flight_light - ground_light) - (flight_dark - ground_dark),
                    levels = design)
fit2 <- eBayes(contrasts.fit(fit, cm))
out <- data.frame(gene = rownames(fit2), check.names = FALSE)
for (k in colnames(cm)) {
  tt <- topTable(fit2, coef = k, number = Inf, sort.by = "none")
  out[[paste0(k, "_logFC")]] <- tt$logFC
  out[[paste0(k, "_t")]] <- tt$t
  out[[paste0(k, "_P")]] <- tt$P.Value
}
write.table(out, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  limma %s, edgeR %s, %d genes, df.prior %.2f\n", packageVersion("limma"),
            packageVersion("edgeR"), nrow(out), fit2$df.prior))
