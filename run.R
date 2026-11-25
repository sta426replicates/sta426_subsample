#!/usr/bin/env Rscript
# Subsample stage: draw n replicates per condition, reproducibly from --seed.
#
# Inputs:  --data.counts (counts.tsv.gz), --data.coldata (coldata.tsv)
# Params:  --n (replicates per condition), --seed
# Outputs (in --output_dir):
#   counts.tsv.gz   same format as the input, restricted to the drawn samples
#   coldata.tsv     same format as the input, restricted to the drawn samples
#   subsample.json  n, seed, RNG settings and the drawn sample ids per condition
#
# All genes are kept, so every downstream DE result covers the same gene set.

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--data.counts", dest = "counts", type = "character"),
  make_option("--data.coldata", dest = "coldata", type = "character"),
  make_option("--n", type = "integer", help = "Replicates per condition"),
  make_option("--seed", type = "integer", help = "Random seed"),
  make_option("--output_dir", type = "character"),
  make_option("--name", type = "character", default = "random_replicates")
)))
for (a in c("counts", "coldata", "n", "seed", "output_dir")) {
  if (is.null(opt[[a]])) stop("Missing required argument: ", a)
}
if (opt$n < 1) stop("--n must be >= 1")
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(gzfile(opt$counts), check.names = FALSE)
coldata <- read.delim(opt$coldata, colClasses = c(sample = "character", condition = "character"))
stopifnot(identical(colnames(counts)[1], "gene"),
          setequal(colnames(counts)[-1], coldata$sample))

conditions <- c("WT", "snf2")
if (!all(coldata$condition %in% conditions)) stop("Unexpected condition labels in coldata")

# Fix the RNG explicitly so a seed gives the same draw on every R version.
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(opt$seed)

drawn <- list()
for (cond in conditions) {
  # Sort first so the draw does not depend on the input column order.
  pool <- sort(coldata$sample[coldata$condition == cond])
  if (opt$n > length(pool)) {
    stop(sprintf("--n %d exceeds the %d %s samples available", opt$n, length(pool), cond))
  }
  drawn[[cond]] <- sort(sample(pool, opt$n))
}
keep <- c(drawn$WT, drawn$snf2)

sub_coldata <- coldata[match(keep, coldata$sample), , drop = FALSE]
sub_counts <- counts[, c("gene", keep), drop = FALSE]

con <- gzfile(file.path(opt$output_dir, "counts.tsv.gz"), "w")
write.table(sub_counts, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
write.table(sub_coldata, file.path(opt$output_dir, "coldata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write_json(list(n = opt$n, seed = opt$seed,
                rng = as.list(RNGkind()),
                samples = drawn),
           file.path(opt$output_dir, "subsample.json"),
           auto_unbox = TRUE, pretty = TRUE)

message(sprintf("n=%d seed=%d: WT {%s} snf2 {%s}", opt$n, opt$seed,
                paste(drawn$WT, collapse = ","), paste(drawn$snf2, collapse = ",")))
