###
##### Simulating DMRs with simDMR function to decide which parameters are the best to detect the greatest number of DMRs

# Set directory
main_dir <- "."
output_dir <- file.path(main_dir, "output")

setwd(main_dir)
require(dmrseq)
require(intervals)
require(BiocParallel)

#parallele processing with MulticoreParam(x), x = number of cores
register(MulticoreParam(6))

### Functions ###

# Get the number of regions in the first set of regions that are sufficiently overlapped by regions in the second set of regions
num_hits <- function(first, second) {
  # Convert regions in first set to intervals
  m <- cbind(start(first), end(first))
  first_intvs <- Intervals(m,  closed = c(TRUE, TRUE), type = "Z")

  # Convert regions in second set to intervals
  m <- cbind(start(second), end(second))
  second_intvs <- Intervals(m,  closed = c(TRUE, TRUE), type = "Z")

  # Determine number of hits
  TH <- 0.2
  n <- 0
  for (i in seq_len(nrow(first_intvs))) {
    intv <- first_intvs[i]
    x <- interval_intersection(intv, second_intvs)
    a <- intv[1, 2] - intv[1, 1] + 1
    b <- if (length(x) == 0) {0} else {x[1, 2] - x[1, 1] + 1}
    if (b / a > TH) {
      n <- n + 1
    }
  }
  return(n)
}

#variables
myminNumRegion <- 3 #default 5
myminInSpan <- 10 #default 30
mycutoff <- 0.01 #default 0.1

# Constants
num_cond <- 2   # number of conditions
num_repl <- c(3, 3)   # number of replicates per condition

# Check number of conditions
if (num_cond != 2) {
  stop("only implemented for num_cond equal to 2")
}

# Check number of replicates
for (i in 1:num_cond) {
  if (num_repl[i] <= 1) {
    stop("only implemented for number of replicates greater than 1 for all conditions")
  }
}

# Read in data
data_dir <- "./data"
counts <- list()   # 2-dim list, first index: condition, second index: replicate
counts[[1]] <- list()
counts[[2]] <- list()

fname <- file.path(data_dir, "MethylKit_F11_EMP_F1_H0_C1.tabular")
counts[[1]][[1]] <- read.csv(fname,   sep = "\t")
fname <- file.path(data_dir, "MethylKit_F11_EMP_F1_H0_C2.tabular")
counts[[1]][[2]]  <- read.csv(fname,  sep = "\t")
fname <- file.path(data_dir, "MethylKit_F11_EMP_F1_H0_C3.tabular")
counts[[1]][[3]]  <- read.csv(fname,  sep = "\t")
fname <- file.path(data_dir, "MethylKit_F11_EMP_F1_H0_NSI1.tabular")
counts[[2]][[1]] <- read.csv(fname,  sep = "\t")
fname <- file.path(data_dir, "MethylKit_F11_EMP_F1_H0_NSI2.tabular")
counts[[2]][[2]] <- read.csv(fname,  sep = "\t")
fname <- file.path(data_dir, "MethylKit_F11_EMP_F1_H0_NSI3.tabular")
counts[[2]][[3]] <- read.csv(fname, sep = "\t")

# Remove columns not needed and rename columns for merging
for (i in 1:num_cond) {
  for (j in 1:num_repl[i]) {
    x <- counts[[i]][[j]]
    y <- x[, c("chr", "base", "strand", "coverage", "freqC")]
    n <- which(colnames(y) == "freqC")
    colnames(y)[n] <- paste("freqC.", i, ".", j, sep = "")
    n <- which(colnames(y) == "coverage")
    colnames(y)[n] <- paste("coverage.", i, ".", j, sep = "")
    counts[[i]][[j]] <- y
  }
}

# Merge count data sets, keeping only positions occurring in all count data sets
l <- list()
for (i in 1:num_cond) {
  l[[i]] <- counts[[i]][[1]]
  for (j in 2:num_repl[i]) {
    l[[i]] <- merge(l[[i]], counts[[i]][[j]], by = c("chr", "base", "strand"))
  }
}
data <- merge(l[[1]], l[[2]], by = c("chr", "base", "strand"))

# Construct BSseq object from count data sets
# ... Construct coverage and methylation frequency
Cov <- c()
M <- c()
for (i in 1:num_cond) {
  for (j in 1:num_repl[i]) {
    col_name <- paste("coverage.", i, ".", j, sep = "")
    Cov <- cbind(Cov, data[, col_name])
    col_name <- paste("freqC.", i, ".", j, sep = "")
    M <- cbind(M, data[, col_name])
  }
}
M <- round(Cov * M / 100)

# ... Contruct BSseq object itself
chr <- data[, "chr"]
pos <- data[, "base"]
sn <- c(paste(1, ".", 1:num_repl[1], sep = ""), paste(2, ".", 1:num_repl[2], sep = ""))
bs <- BSseq(chr = chr, pos = pos, M = M, Cov = Cov, sampleNames = sn)
first_col <- rgb(1, 0.6, 0.6)
second_col <- rgb(0, 0, 1)
ct <- c(rep("first_cond", num_repl[1]), rep("second_cond", num_repl[2]))
rep <- c(paste("replicate", 1:num_repl[1], sep = ""), paste("replicate", 1:num_repl[2], sep = ""))
col <- c(rep(first_col, num_repl[1]), rep(second_col, num_repl[2]))
df <- data.frame(CellType = ct, Rep = rep, col = col)
pData(bs) <- df

# Output number of chromosomes / contigs
v <- as.character(bs@rowRanges@seqnames)
num_chroms <- length(unique(v))
print(paste("number of chromosomes:", num_chroms), quote = FALSE)

# Run dmrseq
bs <- sort(bs)
regions <- dmrseq(bs = bs, cutoff = 0.01, testCovariate = "CellType",
                  minNumRegion = myminNumRegion, minInSpan = myminInSpan, chrsPerChunk = 50)
save(regions, file = file.path(output_dir, "regions_F11_F1_H0_minnumregion3_27.06.19.rda"))

# Plotting
## plotDMRs(bs, regions[4], testCovariate="CellType")


# Export regions as files
# ... Calculate average difference of the methylation degree between the two conditions for each region
meth_perc <- c()
for (i in seq_len(length(regions))) {
  m <- getMeth(bs, regions = regions[i], type = "raw")
  x <- apply(m[[1]], 2, mean)
  a <- mean(x[1:num_repl[1]])
  v <- num_repl[1] + 1:num_repl[2]
  b <- mean(x[v])
  meth_perc <- rbind(meth_perc, c(a, b))
}
diffrc <- meanDiff(bs, dmrs = regions, testCovariate = "CellType")

# ... Prepare data
ranges <- regions@ranges
starts <- ranges@start
widths <- ranges@width
ends <- starts + widths - 1
scaff_names <- as.character(regions@seqnames)
len <- length(starts)

# ... Write to BED file
lines <- c()
for (i in 1:len) {
  l <- paste(scaff_names[i], starts[i], ends[i], sep = "\t")
  lines <- c(lines, l)
}
fc <- file(file.path(output_dir, "regions_F11_F1_H0_minnumregion3_27.06.19.bed"))
writeLines(lines, fc)
close(fc)

write.csv(as.data.frame(regions),  file = file.path(output_dir, "region_F11_F1_H0_minnumregion3_27.06.19.csv"))

# ... Write to custom format file
lines <- c("# chromosome / contig, first position, last position, average methylation degree of the first condition, average methylation degree of the second condition, average difference of the methylation degree between the conditions")
for (i in 1:len) {
  l <- paste(scaff_names[i], starts[i], ends[i], meth_perc[i, 1], meth_perc[i, 2], diffrc[i], sep = "\t")
  lines <- c(lines, l)
}
fc <- file(file.path(output_dir, "regions_percentagediff_F11_F1_H0_minnumregion3_27.06.19.tab"))
writeLines(lines, fc)
close(fc)

write.csv(as.data.frame(fc),
          file = file.path(output_dir, "regions_percentagediff_F11_F1_H0_minnumregion3_27.06.19.csv"))
