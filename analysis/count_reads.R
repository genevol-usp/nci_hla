reads <- commandArgs(TRUE)[1]

n <- length(unique(readLines(reads)))

out <- sub("mapped", "n", reads)

write(n, out)
