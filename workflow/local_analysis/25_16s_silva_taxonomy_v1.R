#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: script.R PROJECT THREADS")
project <- args[[1L]]
threads <- as.integer(args[[2L]])
if (is.na(threads) || threads < 1L) stop("Invalid thread count")

input_dir <- file.path(project, "input", "16S_corrected_v2")
db_dir <- file.path(project, "databases", "silva_138_2_dada2")
out_dir <- file.path(project, "results", "16S_silva1382_v1")
if (dir.exists(out_dir)) stop("Refusing to overwrite existing output directory")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seqtab_file <- file.path(input_dir, "seqtab_nochim.rds")
count_file <- file.path(input_dir, "asv_count_table.tsv")
train_file <- file.path(db_dir, "silva_nr99_v138.2_toGenus_trainset.fa.gz")
species_file <- file.path(db_dir, "silva_v138.2_assignSpecies.fa.gz")
stopifnot(file.exists(seqtab_file), file.exists(count_file),
          file.exists(train_file), file.exists(species_file))

seqtab <- readRDS(seqtab_file)
count_tab <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(is.matrix(seqtab), nrow(seqtab) == 24L, ncol(seqtab) == 9752L,
          nrow(count_tab) == ncol(seqtab), count_tab[[1L]][1L] == "ASV00001")
seqs <- colnames(seqtab)
asv_id <- count_tab[[1L]]

res <- assignTaxonomy(seqs, train_file, minBoot = 80, tryRC = TRUE,
                      multithread = threads, outputBootstraps = TRUE, verbose = TRUE)
tax <- res$tax
boot <- res$boot
tax_species <- addSpecies(tax, species_file, allowMultiple = FALSE, tryRC = TRUE,
                          n = 1e5, verbose = TRUE)

rank_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if (!"Species" %in% colnames(tax_species)) {
  tax_species <- cbind(tax_species, Species = NA_character_)
}
tax_species <- tax_species[, rank_names, drop = FALSE]
boot_out <- matrix(NA_real_, nrow = length(seqs), ncol = length(rank_names),
                   dimnames = list(NULL, rank_names))
boot_out[, colnames(boot)] <- boot

taxonomy <- data.frame(ASV_ID = asv_id, tax_species, check.names = FALSE,
                       stringsAsFactors = FALSE)
bootstrap <- data.frame(ASV_ID = asv_id, boot_out, check.names = FALSE)

kingdom <- ifelse(is.na(taxonomy$Kingdom), "", taxonomy$Kingdom)
family <- ifelse(is.na(taxonomy$Family), "", taxonomy$Family)
genus <- ifelse(is.na(taxonomy$Genus), "", taxonomy$Genus)
is_bacteria_archaea <- kingdom %in% c("Bacteria", "Archaea")
is_chloroplast <- grepl("chloroplast", paste(taxonomy$Class, taxonomy$Order,
                                               family, genus), ignore.case = TRUE)
is_mitochondria <- grepl("mitochond", paste(taxonomy$Class, taxonomy$Order,
                                             family, genus), ignore.case = TRUE)
flags <- data.frame(
  ASV_ID = asv_id,
  is_bacteria_or_archaea = is_bacteria_archaea,
  is_chloroplast = is_chloroplast,
  is_mitochondria = is_mitochondria,
  proposed_keep_for_prokaryote_analysis = is_bacteria_archaea & !is_chloroplast & !is_mitochondria,
  stringsAsFactors = FALSE
)

asv_reads <- rowSums(as.matrix(count_tab[, -1L, drop = FALSE]))
assigned_rank <- function(x) !is.na(x) & nzchar(x)
coverage <- do.call(rbind, lapply(rank_names, function(rank) {
  hit <- assigned_rank(taxonomy[[rank]])
  data.frame(rank = rank, assigned_asvs = sum(hit), asv_fraction = mean(hit),
             assigned_reads = sum(asv_reads[hit]), read_fraction = sum(asv_reads[hit]) / sum(asv_reads))
}))
cleaning <- data.frame(
  category = c("all", "bacteria_or_archaea", "chloroplast", "mitochondria", "proposed_keep"),
  asv_count = c(length(asv_id), sum(is_bacteria_archaea), sum(is_chloroplast),
                sum(is_mitochondria), sum(flags$proposed_keep_for_prokaryote_analysis)),
  read_count = c(sum(asv_reads), sum(asv_reads[is_bacteria_archaea]), sum(asv_reads[is_chloroplast]),
                 sum(asv_reads[is_mitochondria]), sum(asv_reads[flags$proposed_keep_for_prokaryote_analysis]))
)

write.table(taxonomy, file.path(out_dir, "asv_taxonomy_silva1382.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(bootstrap, file.path(out_dir, "asv_taxonomy_bootstrap.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
write.table(flags, file.path(out_dir, "asv_taxonomic_cleaning_flags.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(coverage, file.path(out_dir, "taxonomy_rank_coverage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cleaning, file.path(out_dir, "taxonomy_cleaning_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(tax_species, file.path(out_dir, "taxonomy_silva1382.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
writeLines(c(
  "database=SILVA_138.2_DADA2",
  "taxonomy_method=dada2::assignTaxonomy",
  "minBoot=80", "tryRC=TRUE",
  "species_method=dada2::addSpecies_exact_match",
  "species_allowMultiple=FALSE"
), file.path(out_dir, "taxonomy_parameters.txt"))
cat("TAXONOMY_OK asvs=", length(asv_id), " reads=", sum(asv_reads), "\n", sep = "")


