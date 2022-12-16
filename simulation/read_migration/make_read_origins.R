library(tidyverse)

# Gencode annotations
annotations <- 
    "/raid/genevol/gencode/gencode.v37.primary_assembly.annotation.gtf" |> 
    read_tsv(comment = "#", col_types = "c-cii-c-c",
	     col_names = c("chr", "feature", "start", "end", "strand", "info"))

gene_annot <- annotations |>
    filter(feature == "gene") |>
    mutate(gene_id = str_extract(info, "(?<=gene_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+")) |>
    select(-feature, -info, -strand)

exon_annot <- annotations |>
    filter(feature == "exon") |>
    mutate(transcript_id = str_extract(info, "(?<=transcript_id\\s\")[^\"]+"),
	   gene_name = str_extract(info, "(?<=gene_name\\s\")[^\"]+")) |>
    select(-feature, -info) 

write_tsv(gene_annot, "gene_coords.tsv")

# Simulated reads for each individual
samples <- read_lines("../samples.txt")

read_ids <- "/scratch/vitor/temp/%s_readids.txt" |>
    sprintf(samples) |>
    setNames(samples) |>
    map_dfr(function(x) read_lines(x) |>
	    tibble(id = _) |>
	    separate(id, c("id", "info"), sep = " ") |>
	    separate(info, c("read", "mate"), sep = ";") |>
	    mutate(id = sub("^@", "", id),
		   read = sub("mate1:", "", read),
		   mate = sub("mate2:", "", mate),
		   mate = sub("/1$", "", mate)) |>
	    separate(read, c("read1", "read2"), sep = "-", convert = TRUE) |>
	    separate(mate, c("mate1", "mate2"), sep = "-", convert = TRUE),
	    .id = "sampleid") |>
    mutate(transcript_id = sub("^read\\d+_((ENST[0-9.]+)(_.+)?)$", "\\1", id))
    
# Find origin position for REF
read_ids_ref <- filter(read_ids, !grepl("_[ABC]\\*", id)) |>
    left_join(distinct(exon_annot, transcript_id, strand)) |>
    mutate(read = ifelse(strand == "+", read1, read2),
	   mate = ifelse(strand == "+", mate1, mate2)) |>
    select(sampleid, readid = id, transcript_id, read, mate)

ref_transcripts_positions <- exon_annot |>
    filter(transcript_id %in% read_ids_ref$transcript_id) |>
    mutate(positions = map2(start, end, ~.x:.y)) |>
    select(chr, transcript_id, gene_name, strand, positions) |>
    unnest(cols = c(positions))

ref_transcripts_positions_plus <- ref_transcripts_positions |>
    filter(strand == "+") |>
    group_split(transcript_id) |>
    map_dfr(~arrange(., positions) |> rowid_to_column("ix"))

ref_transcripts_positions_minus <- ref_transcripts_positions |>
    filter(strand == "-") |>
    group_split(transcript_id) |>
    map_dfr(~arrange(., desc(positions)) |> rowid_to_column("ix"))

ref_positions_df <- bind_rows(ref_transcripts_positions_plus, ref_transcripts_positions_minus) |>
    select(-strand)

ref_reads_origins <- inner_join(read_ids_ref, ref_positions_df,
	   join_by(transcript_id, read == ix)) |>
    select(sampleid, readid, chr, gene_name, origin_read = positions)

ref_mates_origins <- inner_join(read_ids_ref, ref_positions_df,
	   join_by(transcript_id, mate == ix)) |>
    select(sampleid, readid, chr, gene_name, origin_mate = positions)

ref_origins <- left_join(ref_reads_origins, ref_mates_origins) |>
    select(sampleid, readid, chr, gene_name, origin_read, origin_mate)

# Find origin position for HLA
hladb <- read_tsv("../../indices/personalize_transcripts/position_personalized_transcripts.tsv")

read_ids_hla <- filter(read_ids, grepl("_[ABC]\\*", id)) |>
    mutate(gene_name = sub("^read\\d+_ENST[0-9.]+_([ABC])\\*.+$", "HLA-\\1", id),
	   read = ifelse(gene_name == "HLA-A", read1, read2),
	   mate = ifelse(gene_name == "HLA-A", mate1, mate2)) |>
    select(sampleid, readid = id, gene_name, transcript_id, read, mate)

hla_reads_origins <- 
    inner_join(read_ids_hla, hladb, 
	   join_by(transcript_id == tx, read == i)) %>%
    select(sampleid, readid, gene_name, origin_read = pos)

hla_mates_origins <- 
    inner_join(read_ids_hla, hladb, 
	   join_by(transcript_id == tx, mate == i)) %>%
    select(sampleid, readid, gene_name, origin_mate = pos)

hla_origins <- left_join(hla_reads_origin, hla_mates_origin) |>
    mutate(chr = "chr6") |>
    select(sampleid, readid, chr, gene_name, origin_read, origin_mate)


origins_out <- bind_rows(ref_origins, hla_origins) |>
    arrange(sampleid, readid)

write_tsv(origins_out, "read_origins.tsv")
