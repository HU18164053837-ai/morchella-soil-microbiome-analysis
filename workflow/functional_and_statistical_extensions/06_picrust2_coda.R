args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript 06_picrust2_coda.R <pathways.tsv.gz> <weighted_nsti.tsv.gz> <outdir>")
path_file <- args[1]; nsti_file <- args[2]; outdir <- args[3]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tab <- read.delim(gzfile(path_file), check.names = FALSE, row.names = 1)
samples <- colnames(tab)
ma <- samples[grepl("^[MA][234]-[123]$", samples)]
pairs <- unlist(lapply(c(2,3,4), function(g) paste0(g, "-", 1:3)))
ordered <- c(paste0("M", pairs), paste0("A", pairs))
if (!all(ordered %in% ma)) stop("Expected nine matched M/A pairs were not found")
x <- as.matrix(tab[, ordered, drop = FALSE])
storage.mode(x) <- "numeric"

meta <- data.frame(sample = ordered,
                   stage = factor(rep(c("M", "A"), each = 9), levels = c("M", "A")),
                   pair_id = factor(rep(pairs, 2)),
                   greenhouse = factor(rep(paste0("GH", substr(pairs,1,1)), 2)),
                   row.names = ordered)

nsti <- read.delim(gzfile(nsti_file), check.names = FALSE)
names(nsti)[1:2] <- c("sample", "weighted_NSTI")
meta$weighted_NSTI <- nsti$weighted_NSTI[match(rownames(meta), nsti$sample)]
write.table(meta, file.path(outdir, "MA_sample_metadata_with_NSTI.tsv"), sep = "\t", quote = FALSE, col.names = NA)

bh <- function(p) p.adjust(p, method = "BH")
exact_signflip <- function(d) {
  obs <- mean(d)
  signs <- as.matrix(expand.grid(rep(list(c(-1,1)), length(d))))
  perm <- as.vector(signs %*% d / length(d))
  c(effect = obs, p = (1 + sum(abs(perm) >= abs(obs) - 1e-15))/(1+nrow(signs)), positive = sum(d > 0))
}

# Design-aware exploratory CLR analysis. A small relative-abundance pseudocount is
# used only after sample closure; sensitivity is evaluated separately below.
rel <- sweep(x, 2, colSums(x), "/")
clr_one <- function(pc) {
  z <- log(rel + pc)
  z <- sweep(z, 2, colMeans(z), "-")
  out <- t(vapply(seq_len(nrow(z)), function(i) {
    d <- z[i, paste0("A", pairs)] - z[i, paste0("M", pairs)]
    e <- exact_signflip(d)
    loo <- sapply(c("2","3","4"), function(g) mean(d[substr(pairs,1,1) != g]))
    c(mean_delta = unname(e["effect"]),
      exact_p = unname(e["p"]),
      pairs_positive = unname(e["positive"]),
      loo_GH2 = unname(loo[1]),
      loo_GH3 = unname(loo[2]),
      loo_GH4 = unname(loo[3]))
  }, numeric(6)))
  out <- data.frame(pathway = rownames(z), pseudocount = pc, out, check.names = FALSE)
  out$BH_q <- bh(out$exact_p)
  out
}
clr_all <- do.call(rbind, lapply(c(1e-8, 1e-7, 1e-6), clr_one))
stopifnot(all(c("pathway", "pseudocount", "mean_delta", "exact_p", "BH_q",
                "pairs_positive", "loo_GH2", "loo_GH3", "loo_GH4") %in%
              names(clr_all)))
write.table(clr_all, file.path(outdir, "MetaCyc_CLR_exact_paired_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ALDEx2 requires count-like non-negative inputs; PICRUSt2 inferred abundances are
# rounded only for its Dirichlet Monte Carlo layer. Raw predictions remain archived.
suppressPackageStartupMessages(library(ALDEx2))
ald_counts <- round(x)
ald_groups <- as.character(meta$stage)
stopifnot(
  is.matrix(ald_counts),
  is.numeric(ald_counts),
  length(ald_groups) == ncol(ald_counts),
  is.character(ald_groups),
  all(ald_groups %in% c("M", "A"))
)
ald <- aldex.clr(ald_counts, ald_groups, mc.samples = 128, denom = "all", verbose = FALSE)
ald_t <- aldex.ttest(ald, paired.test = TRUE, verbose = FALSE)
ald_e <- aldex.effect(ald, verbose = FALSE)
ald_out <- data.frame(pathway = rownames(ald_t), ald_t, ald_e[rownames(ald_t), , drop = FALSE], check.names = FALSE)
write.table(ald_out, file.path(outdir, "MetaCyc_ALDEx2_paired.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Dependency-free compositional consensus and NSTI audit. ANCOMBC 2.8.0 requires
# phyloseq/TreeSummarizedExperiment, which are absent from the locked offline
# environment. We therefore triangulate paired ALDEx2 with exact paired CLR,
# pseudocount sensitivity and leave-one-greenhouse-out direction stability.
ald_q_col <- if ("we.eBH" %in% names(ald_out)) "we.eBH" else grep("BH", names(ald_out), value=TRUE)[1]
ald_eff_col <- if ("effect" %in% names(ald_out)) "effect" else grep("effect", names(ald_out), value=TRUE)[1]
clr_ref <- subset(clr_all, pseudocount == 1e-7)
cons <- merge(
  data.frame(pathway=ald_out$pathway,
             aldex_q=ald_out[[ald_q_col]],
             aldex_effect_raw_M_minus_A=ald_out[[ald_eff_col]],
             aldex_effect_A_minus_M=-ald_out[[ald_eff_col]]),
  clr_ref[,c("pathway","mean_delta","exact_p","BH_q","pairs_positive",
            "loo_GH2","loo_GH3","loo_GH4")],
  by="pathway", all=TRUE
)
names(cons)[names(cons)=="BH_q"] <- "clr_exact_BH_q"
cons$aldex_clr_direction_agree <- sign(cons$aldex_effect_A_minus_M) == sign(cons$mean_delta)
cons$loo_direction_stable <- with(cons, sign(mean_delta)==sign(loo_GH2) & sign(mean_delta)==sign(loo_GH3) & sign(mean_delta)==sign(loo_GH4))
cons$consensus_q10 <- with(cons,
  !is.na(aldex_q) & !is.na(clr_exact_BH_q) &
  aldex_q < 0.10 & clr_exact_BH_q < 0.10 &
  aldex_clr_direction_agree & loo_direction_stable
)

delta_nsti <- meta[paste0("A",pairs),"weighted_NSTI"] - meta[paste0("M",pairs),"weighted_NSTI"]
nsti_test <- exact_signflip(delta_nsti)
nsti_report <- data.frame(mean_M=mean(meta[paste0("M",pairs),"weighted_NSTI"]),
                          mean_A=mean(meta[paste0("A",pairs),"weighted_NSTI"]),
                          mean_delta=nsti_test["effect"], pairs_positive=nsti_test["positive"],
                          exact_p=nsti_test["p"])
write.table(nsti_report, file.path(outdir, "NSTI_paired_stage_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Association between each predicted pathway change and NSTI change is a warning flag,
# not a correction or proof of bias.
pc <- 1e-7; z <- log(rel+pc); z <- sweep(z,2,colMeans(z),"-")
nsti_cor <- t(vapply(seq_len(nrow(z)), function(i) {
  d <- z[i,paste0("A",pairs)]-z[i,paste0("M",pairs)]
  ct <- suppressWarnings(cor.test(d, delta_nsti, method="spearman", exact=FALSE))
  c(rho=unname(ct$estimate), p=ct$p.value)
}, numeric(2)))
nsti_cor <- data.frame(pathway=rownames(z), nsti_delta_rho=nsti_cor[,1], nsti_delta_p=nsti_cor[,2])
nsti_cor$nsti_delta_q <- bh(nsti_cor$nsti_delta_p)
cons <- merge(cons, nsti_cor, by="pathway", all.x=TRUE)
cons$nsti_sensitive_q10 <- !is.na(cons$nsti_delta_q) & cons$nsti_delta_q < 0.10
cons$paper_candidate <- cons$consensus_q10 & cons$loo_direction_stable & !cons$nsti_sensitive_q10
cons <- cons[order(!cons$paper_candidate,
                   pmax(cons$aldex_q, cons$clr_exact_BH_q, na.rm=TRUE)), ]
write.table(cons, file.path(outdir, "MetaCyc_consensus_and_NSTI_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

sink(file.path(outdir, "PICRUSt2_CoDA_REPORT.txt"))
cat("PICRUSt2 predicted-function CoDA analysis complete\n")
cat("Pathways:", nrow(x), " Samples:", ncol(x), " M/A pairs: 9\n")
cat("ALDEx2 q<0.10:", sum(cons$aldex_q < 0.10, na.rm=TRUE), "\n")
cat("Exact paired CLR q<0.10:", sum(cons$clr_exact_BH_q < 0.10, na.rm=TRUE), "\n")
cat("ALDEx2/CLR robust consensus q<0.10:", sum(cons$consensus_q10, na.rm=TRUE), "\n")
cat("Final candidates after leave-one-greenhouse and NSTI audit:", sum(cons$paper_candidate, na.rm=TRUE), "\n")
cat("Boundary: PICRUSt2 is predicted genomic potential derived from 16S taxonomy, not measured genes or activity.\n")
sink()

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
cat("MORCHELLA_PICRUST2_CODA_COMPLETE\n")

