#!/usr/bin/env Rscript

.libPaths(c("C:/path/to/R/library", .libPaths()))
suppressPackageStartupMessages(library(vegan))

root <- "C:/path/to/morchella_microbiome/analysis"
out <- file.path(root, "stageA_reviewer_enhancement_20260823_v5")
if (dir.exists(out)) stop("Refusing to overwrite existing output: ", out)
dir.create(out, recursive = TRUE)
dir.create(file.path(out, "greenhouse_level"))
dir.create(file.path(out, "cross_domain"))
dir.create(file.path(out, "picrust2"))
dir.create(file.path(out, "funguild"))
dir.create(file.path(out, "candidate_audit"))

write_tsv <- function(x, path) write.table(x, path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
bh <- function(p) p.adjust(p, "BH")
exact_signflip <- function(d) {
  d <- as.numeric(d); n <- length(d); obs <- mean(d)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  perm <- as.numeric(signs %*% d / n)
  c(mean_delta=obs, exact_p=(1 + sum(abs(perm) >= abs(obs) - 1e-15))/(1 + length(perm)),
    positive=sum(d > 0), negative=sum(d < 0))
}
sample_order <- unlist(lapply(c(2,3,4), function(g) paste0(g, "-", 1:3)))
M_ids <- paste0("M", sample_order); A_ids <- paste0("A", sample_order)
gh <- paste0("GH", substr(sample_order, 1, 1))

# ---------- 1. Greenhouse-level summaries ----------
paired <- read.delim(file.path(root, "statistical_enhancement_20260823", "paired_effect_input_long.tsv"), check.names=FALSE)
gh_summary <- aggregate(delta_A_minus_M ~ family + marker + feature + unit + greenhouse,
                        data=paired, FUN=mean)
names(gh_summary)[names(gh_summary)=="delta_A_minus_M"] <- "greenhouse_mean_A_minus_M"
direction_summary <- do.call(rbind, lapply(split(gh_summary, interaction(gh_summary$family, gh_summary$marker, gh_summary$feature, drop=TRUE)), function(z) {
  data.frame(family=z$family[1], marker=z$marker[1], feature=z$feature[1], unit=z$unit[1],
             greenhouses_positive=sum(z$greenhouse_mean_A_minus_M > 0),
             greenhouses_negative=sum(z$greenhouse_mean_A_minus_M < 0),
             greenhouse_direction_consistent=all(z$greenhouse_mean_A_minus_M > 0) || all(z$greenhouse_mean_A_minus_M < 0),
             min_greenhouse_effect=min(z$greenhouse_mean_A_minus_M), max_greenhouse_effect=max(z$greenhouse_mean_A_minus_M))
}))
write_tsv(gh_summary, file.path(out, "greenhouse_level", "all_existing_effects_by_greenhouse.tsv"))
write_tsv(direction_summary, file.path(out, "greenhouse_level", "all_existing_effects_greenhouse_direction_summary.tsv"))

# Community within-pair Bray-Curtis shifts.
read_asv <- function(path) {
  z <- read.delim(path, check.names=FALSE)
  if ("sample_id" %in% names(z)) {
    m <- as.matrix(z[,setdiff(names(z),"sample_id"),drop=FALSE]); storage.mode(m)<-"numeric"
    rownames(m)<-z$sample_id
  } else {
    id <- if ("ASV_ID" %in% names(z)) "ASV_ID" else names(z)[1]
    m <- t(as.matrix(z[, setdiff(names(z), id), drop=FALSE])); storage.mode(m) <- "numeric"
    colnames(m) <- z[[id]]
  }
  m
}
cnt16 <- read_asv(file.path(root, "downstream", "16S_v1_silva_diversity", "objects", "16S_rarefied_counts_depth1500.tsv"))
cntits <- read_asv(file.path(root, "frozen", "ITS_STAGE1_FREEZE_2026-08-21_v1", "evidence", "analysis", "downstream", "ITS_v2_taxonomy", "tables", "fungal_asv_count_table.tsv"))
set.seed(20260823L)
cntits <- rrarefy(cntits, 16000)
d16 <- as.matrix(vegdist(cnt16, "bray")); dits <- as.matrix(vegdist(cntits, "bray"))
pair_shift <- rbind(
  data.frame(marker="16S", pair_id=sample_order, greenhouse=gh, within_pair_bray=mapply(function(m,a)d16[m,a], M_ids,A_ids)),
  data.frame(marker="ITS", pair_id=sample_order, greenhouse=gh, within_pair_bray=mapply(function(m,a)dits[m,a], M_ids,A_ids)))
write_tsv(pair_shift, file.path(out, "greenhouse_level", "community_within_pair_bray_shift.tsv"))
pair_shift_gh <- aggregate(within_pair_bray ~ marker + greenhouse, pair_shift, mean)
write_tsv(pair_shift_gh, file.path(out, "greenhouse_level", "community_within_pair_bray_by_greenhouse.tsv"))

# ---------- 2. Cross-domain structure sensitivity ----------
ma_ids <- c(M_ids, A_ids)
meta <- data.frame(sample=ma_ids, stage=substr(ma_ids,1,1), pair=sub("^[MA]", "", ma_ids),
                   greenhouse=paste0("GH", substr(sub("^[MA]", "", ma_ids),1,1)))
mantel_one <- function(ids, label, blocks=NULL, partial_design=FALSE) {
  x <- as.dist(d16[ids,ids]); y <- as.dist(dits[ids,ids])
  if (is.null(blocks)) {
    fit <- mantel(x,y,method="spearman",permutations=9999)
  } else {
    ctrl <- permute::how(nperm=9999, blocks=factor(blocks))
    fit <- mantel(x,y,method="spearman",permutations=ctrl)
  }
  data.frame(analysis=label, n_samples=length(ids), mantel_r=unname(fit$statistic), p_value=fit$signif,
             permutation_structure=ifelse(is.null(blocks),"unrestricted","within_greenhouse"))
}
cross_tests <- rbind(
  mantel_one(ma_ids, "all_MA_reference_unrestricted"),
  mantel_one(M_ids, "within_M_stage", gh),
  mantel_one(A_ids, "within_A_stage", gh)
)
# Paired community-change magnitude correlation plus leave-one-greenhouse-out.
wide_shift <- merge(pair_shift[pair_shift$marker=="16S",], pair_shift[pair_shift$marker=="ITS",], by=c("pair_id","greenhouse"), suffixes=c("_16S","_ITS"))
rho_all <- cor.test(wide_shift$within_pair_bray_16S, wide_shift$within_pair_bray_ITS, method="spearman", exact=FALSE)
rho_loo <- do.call(rbind, lapply(unique(wide_shift$greenhouse), function(g) {
  z <- wide_shift[wide_shift$greenhouse != g,]
  data.frame(analysis=paste0("paired_change_magnitude_leave_out_",g), n_pairs=nrow(z),
             spearman_rho=cor(z$within_pair_bray_16S,z$within_pair_bray_ITS,method="spearman"))
}))
rho_main <- data.frame(analysis="paired_change_magnitude_all", n_pairs=nrow(wide_shift), spearman_rho=unname(rho_all$estimate), p_value=rho_all$p.value)
write_tsv(cross_tests, file.path(out, "cross_domain", "stage_stratified_mantel.tsv"))
write_tsv(wide_shift, file.path(out, "cross_domain", "paired_cross_domain_change_magnitudes.tsv"))
write_tsv(rbind(transform(rho_main, p_value=p_value), transform(rho_loo, p_value=NA_real_)), file.path(out, "cross_domain", "paired_change_magnitude_correlations.tsv"))

# ---------- 3. PICRUSt2 modules without pathway preselection ----------
ab <- read.delim(gzfile(file.path(root, "extension_hpc_20260822", "results", "picrust2_coda_v1", "MetaCyc_path_abundance_raw.tsv.gz")), check.names=FALSE, row.names=1)
annot <- read.delim(file.path(root, "extension_hpc_20260822", "results", "picrust2_core_v1", "MetaCyc_candidates_annotated_ranked.tsv"), check.names=FALSE)
keep_samples <- c(M_ids,A_ids)
x <- as.matrix(ab[,keep_samples,drop=FALSE]); storage.mode(x)<-"numeric"
rel <- sweep(x,2,colSums(x),"/"); clr <- log(rel+1e-7); clr <- sweep(clr,2,colMeans(clr),"-")
annot <- annot[!duplicated(annot$pathway) & annot$pathway %in% rownames(clr),]
mods <- sort(unique(annot$module[annot$module != "Other predicted metabolism"]))
module_scores <- do.call(rbind,lapply(mods,function(mod){
  ids <- annot$pathway[annot$module==mod]
  sc <- apply(clr[ids,,drop=FALSE],2,median)
  data.frame(module=mod,n_all_pathways=length(ids),sample=names(sc),score=as.numeric(sc))
}))
module_delta <- do.call(rbind,lapply(split(module_scores,module_scores$module),function(z){
  data.frame(module=z$module[1],n_all_pathways=z$n_all_pathways[1],pair_id=sample_order,greenhouse=gh,
             delta_A_minus_M=z$score[match(A_ids,z$sample)]-z$score[match(M_ids,z$sample)])
}))
module_stats <- do.call(rbind,lapply(split(module_delta,module_delta$module),function(z){
  ex <- exact_signflip(z$delta_A_minus_M)
  gm <- aggregate(delta_A_minus_M~greenhouse,z,mean)
  data.frame(module=z$module[1],n_all_pathways=z$n_all_pathways[1],mean_delta=ex["mean_delta"],median_delta=median(z$delta_A_minus_M),
             positive_pairs=ex["positive"],exact_p=ex["exact_p"],greenhouses_positive=sum(gm$delta_A_minus_M>0),
             greenhouses_negative=sum(gm$delta_A_minus_M<0),greenhouse_direction_consistent=all(gm$delta_A_minus_M>0)||all(gm$delta_A_minus_M<0))
}))
module_stats$BH_q <- bh(module_stats$exact_p)
module_gh <- aggregate(delta_A_minus_M~module+n_all_pathways+greenhouse,module_delta,mean)
nsti <- read.delim(file.path(root,"extension_hpc_20260822","results","picrust2_coda_v1","MA_sample_metadata_with_NSTI.tsv"),check.names=FALSE)
nsti_delta <- data.frame(pair_id=sample_order,greenhouse=gh,M=nsti$weighted_NSTI[match(M_ids,nsti$sample)],A=nsti$weighted_NSTI[match(A_ids,nsti$sample)])
nsti_delta$delta_A_minus_M <- nsti_delta$A-nsti_delta$M
nsti_gh <- aggregate(cbind(M,A,delta_A_minus_M)~greenhouse,nsti_delta,mean)
high_pair <- nsti_delta$pair_id[which.max(nsti_delta$A)]
module_nsti_sens <- do.call(rbind,lapply(split(module_delta,module_delta$module),function(z){
  keep <- z$pair_id != high_pair
  data.frame(module=z$module[1],full_mean=mean(z$delta_A_minus_M),mean_excluding_highest_A_NSTI_pair=mean(z$delta_A_minus_M[keep]),
             direction_preserved=sign(mean(z$delta_A_minus_M))==sign(mean(z$delta_A_minus_M[keep])),excluded_pair=high_pair,
             spearman_delta_vs_NSTI=cor(z$delta_A_minus_M,nsti_delta$delta_A_minus_M[match(z$pair_id,nsti_delta$pair_id)],method="spearman"))
}))
write_tsv(module_scores,file.path(out,"picrust2","all_pathway_module_CLR_scores.tsv"))
write_tsv(module_delta,file.path(out,"picrust2","all_pathway_module_paired_deltas.tsv"))
write_tsv(module_stats,file.path(out,"picrust2","all_pathway_module_statistics.tsv"))
write_tsv(module_gh,file.path(out,"picrust2","all_pathway_module_greenhouse_means.tsv"))
write_tsv(nsti_delta,file.path(out,"picrust2","NSTI_pair_level.tsv"))
write_tsv(nsti_gh,file.path(out,"picrust2","NSTI_greenhouse_level.tsv"))
write_tsv(module_nsti_sens,file.path(out,"picrust2","module_NSTI_sensitivity.tsv"))

# ---------- 4. FUNGuild coverage and dominant-taxon sensitivity ----------
its_all <- read.delim(file.path(root,"extension_hpc_20260822","input","ITS_asv_count_table_remote_used.tsv"),check.names=FALSE)
names(its_all) <- sub("_ITS_trimmed\\.fastq\\.gz$", "", names(its_all))
ann <- read.delim(file.path(root,"hpc_results","extension_analysis_20260822","funguild_v3","ITS_ASV_FUNGuild_offline.tsv"),check.names=FALSE,quote="")
idcol <- if("ASV_ID" %in% names(its_all)) "ASV_ID" else names(its_all)[1]
ann <- ann[match(its_all[[idcol]],ann$ASV_ID),]
countmat <- as.matrix(its_all[,keep_samples]); storage.mode(countmat)<-"numeric"; rownames(countmat)<-its_all[[idcol]]
split_mode <- function(s) { if(is.na(s)||!nzchar(s)) character() else trimws(unlist(strsplit(s,"-",fixed=TRUE))) }
allowed_sets <- list(strict=c("Probable","Highly Probable"),inclusive=c("Possible","Probable","Highly Probable"))
calc_fg <- function(allowed, excluded_genus=character(), excluded_asv=character()) {
  ok <- ann$confidence_ranking %in% allowed & !(ann$Genus %in% excluded_genus) & !(ann$ASV_ID %in% excluded_asv)
  modes <- lapply(ann$trophic_mode,split_mode)
  ok <- ok & lengths(modes)>0
  annotated_total <- colSums(countmat[ok,,drop=FALSE])
  sap <- rep(0,ncol(countmat)); names(sap)<-colnames(countmat)
  for(i in which(ok)) if("Saprotroph" %in% modes[[i]]) sap <- sap + countmat[i,]/length(modes[[i]])
  all_total <- colSums(countmat)
  vals <- list(all_reads=sap/all_total,annotated_reads=sap/annotated_total)
  do.call(rbind,lapply(names(vals),function(den){
    v<-vals[[den]]; d<-v[A_ids]-v[M_ids]; ex<-exact_signflip(d); gm<-tapply(d,gh,mean)
    data.frame(denominator=den,mean_M=mean(v[M_ids]),mean_A=mean(v[A_ids]),mean_delta=ex["mean_delta"],positive_pairs=ex["positive"],exact_p=ex["exact_p"],
               GH2=gm["GH2"],GH3=gm["GH3"],GH4=gm["GH4"],greenhouse_direction_consistent=all(gm>0)||all(gm<0),
               mean_annotated_read_coverage_M=mean(annotated_total[M_ids]/all_total[M_ids]),mean_annotated_read_coverage_A=mean(annotated_total[A_ids]/all_total[A_ids]))
  }))
}
# Identify the largest strict saprotroph contributor after cultivation.
modes_all <- lapply(ann$trophic_mode,split_mode)
strict_ok <- ann$confidence_ranking %in% allowed_sets$strict & lengths(modes_all)>0
sap_contrib_A <- sapply(seq_len(nrow(ann)),function(i) if(strict_ok[i] && "Saprotroph" %in% modes_all[[i]]) sum(countmat[i,A_ids])/length(modes_all[[i]]) else 0)
top_asv <- ann$ASV_ID[which.max(sap_contrib_A)]
scenarios <- list(none=list(g=character(),a=character()),remove_Mortierella=list(g="Mortierella",a=character()),
                  remove_Botryotrichum=list(g="Botryotrichum",a=character()),remove_both_genera=list(g=c("Mortierella","Botryotrichum"),a=character()),
                  remove_top_saprotroph_ASV=list(g=character(),a=top_asv))
fg_sens <- do.call(rbind,lapply(names(allowed_sets),function(cs) do.call(rbind,lapply(names(scenarios),function(sn){
  z<-calc_fg(allowed_sets[[cs]],scenarios[[sn]]$g,scenarios[[sn]]$a); z$confidence_set<-cs; z$scenario<-sn; z$excluded_asv<-ifelse(sn=="remove_top_saprotroph_ASV",top_asv,""); z
}))))
fg_sens <- fg_sens[,c("confidence_set","scenario","excluded_asv",setdiff(names(fg_sens),c("confidence_set","scenario","excluded_asv")))]
coverage <- do.call(rbind,lapply(names(allowed_sets),function(cs){
  ok<-ann$confidence_ranking %in% allowed_sets[[cs]] & lengths(modes_all)>0
  total<-colSums(countmat); ac<-colSums(countmat[ok,,drop=FALSE])
  data.frame(confidence_set=cs,sample_id=keep_samples,stage=substr(keep_samples,1,1),pair_id=sub("^[MA]","",keep_samples),greenhouse=paste0("GH",substr(sub("^[MA]","",keep_samples),1,1)),total_reads=total,annotated_reads=ac,annotated_read_fraction=ac/total)
}))
write_tsv(coverage,file.path(out,"funguild","annotation_coverage_by_sample.tsv"))
write_tsv(fg_sens,file.path(out,"funguild","saprotroph_dominant_taxon_sensitivity.tsv"))
write_tsv(data.frame(top_strict_saprotroph_ASV=top_asv,genus=ann$Genus[match(top_asv,ann$ASV_ID)],A_saprotroph_fractional_reads=max(sap_contrib_A)),file.path(out,"funguild","top_saprotroph_ASV.tsv"))

# ---------- 5. Candidate freezing and external grading audit ----------
freeze <- read.delim(file.path(root,"extension_hpc_20260822","input","frozen_candidates.tsv"),check.names=FALSE)
audit <- read.delim(file.path(root,"core_genera_audit_20260822","eight_genera_audit_summary.tsv"),check.names=FALSE)
external <- read.delim(file.path(root,"public_validation","results","final_cross_project_evidence_v3","final_cross_project_evidence_matrix.tsv"),check.names=FALSE)
audit8 <- Reduce(function(x,y) merge(x,y,by=c("marker","genus"),all=TRUE,suffixes=c("","_dup")),list(freeze,audit,external))
audit8$freeze_timing <- "frozen_before_external_dataset_inspection"
audit8$local_inference_status <- "discovery_selected; local p/q values are post-selection descriptive evidence"
audit8$external_inference_status <- "directional_stress_test_not_validation_cohort"
write_tsv(audit8,file.path(out,"candidate_audit","Table_S2_frozen_candidates_full_audit.tsv"))
scr16 <- read.delim(file.path(root,"downstream","16S_v2_compositional_soil","tables","MA_genus_compositional_screening.tsv"),check.names=FALSE)
scr16$marker<-"16S"; scr16$frozen_candidate<-scr16$genus %in% freeze$genus[freeze$marker=="16S"]
scrits <- read.delim(file.path(root,"downstream","ITS_v4_genus_screening","tables","MA_genus_screening.tsv"),check.names=FALSE)
scrits$marker<-"ITS"; scrits$frozen_candidate<-scrits$genus %in% freeze$genus[freeze$marker=="ITS"]
write_tsv(scr16,file.path(out,"candidate_audit","initial_candidate_pool_16S.tsv"))
write_tsv(scrits,file.path(out,"candidate_audit","initial_candidate_pool_ITS.tsv"))
ext_effects <- read.delim(file.path(root,"public_validation","results","final_cross_project_evidence_v3","final_external_effect_sizes.tsv"),check.names=FALSE)
write_tsv(ext_effects,file.path(out,"candidate_audit","external_contrast_level_evidence.tsv"))

# ---------- Human-readable report and manifest ----------
sink(file.path(out,"STAGE_A_RESULTS_SUMMARY.txt"))
cat("STAGE A REVIEWER ENHANCEMENT COMPLETE\n")
cat("Created:",format(Sys.time()),"\n\n")
cat("Cross-domain stage-stratified Mantel results\n"); print(cross_tests,row.names=FALSE)
cat("\nPaired change-magnitude correlation\n"); print(rho_main,row.names=FALSE)
cat("\nPICRUSt2 all-pathway module statistics\n"); print(module_stats,row.names=FALSE)
cat("\nNSTI greenhouse means\n"); print(nsti_gh,row.names=FALSE)
cat("\nFUNGuild saprotroph sensitivity\n"); print(fg_sens,row.names=FALSE)
cat("\nTop saprotroph ASV:",top_asv,"; genus",ann$Genus[match(top_asv,ann$ASV_ID)],"\n")
sink()
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
files <- list.files(out,recursive=TRUE,full.names=TRUE)
manifest <- data.frame(file=sub(paste0("^",gsub("([\\W])","\\\\\\1",out),"[/\\\\]?"),"",files),bytes=file.info(files)$size)
write_tsv(manifest,file.path(out,"OUTPUT_MANIFEST.tsv"))
cat("STAGE_A_COMPLETE\n",out,"\n")

