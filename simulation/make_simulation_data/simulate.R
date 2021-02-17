library(polyester)

col_i <- as.integer(commandArgs(TRUE)[1])

phenotypes <- as.matrix(readr::read_tsv("./phenotypes.tsv")[, col_i])

sample_i <- colnames(phenotypes)
outdir <- file.path("/scratch/vitor/simulation", sample_i) 

dir.create(outdir)

simulate_experiment_countmat(fasta = "./simulation_index.fa",
			     readmat = phenotypes,
			     readlen = 126,
			     fraglen = 300,
			     outdir = outdir)

system(sprintf("./fa2fq.sh %s && rm %s/*fasta", outdir, outdir))
