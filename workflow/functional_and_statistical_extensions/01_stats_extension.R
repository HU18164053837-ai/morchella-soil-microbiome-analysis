options(stringsAsFactors=FALSE)
args<-commandArgs(TRUE); base<-args[1]; out<-args[2]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
rd<-function(x) read.delim(x,check.names=FALSE,quote="",comment.char="")
meta<-rd(file.path(base,"metadata/sample_metadata.tsv")); if(!"sample_id"%in%names(meta)) meta$sample_id<-meta$biological_sample_id
soil<-rd(file.path(base,"extension_analysis_20260822/input/soil_metrics_normalized_v2.tsv"))
cand<-rd(file.path(base,"extension_analysis_20260822/input/frozen_candidates.tsv"))
analyze_marker<-function(marker,countfile,taxfile){
  ct<-rd(countfile); tx<-rd(taxfile); rownames(ct)<-ct[[1]]; ct<-as.matrix(ct[,-1]); storage.mode(ct)<-"numeric"
  colnames(ct)<-sub("_(16S|ITS)_trimmed\\.fastq\\.gz$","",colnames(ct))
  rownames(tx)<-tx[[1]]; tx<-tx[rownames(ct),,drop=FALSE]
  sm<-meta[meta$marker==marker & meta$cohort=="cultivation_paired",]; sm<-sm[sm$sample_id %in% colnames(ct),]
  sm<-sm[match(c(paste0("M",rep(2:4,each=3),"-",1:3),paste0("A",rep(2:4,each=3),"-",1:3)),sm$sample_id),]; sm<-sm[!is.na(sm$sample_id),]
  x<-t(ct[,sm$sample_id,drop=FALSE]); rel<-x/rowSums(x)
  bray<-function(z){n<-nrow(z); d<-matrix(0,n,n); for(i in 1:(n-1)) for(j in (i+1):n) d[i,j]<-d[j,i]<-sum(abs(z[i,]-z[j,]))/sum(z[i,]+z[j,]); d}
  D<-bray(rel); stage<-ifelse(substr(sm$sample_id,1,1)=="A",1,0); pair<-sm$pair_id
  pseudoF<-function(g,dmat=D){gm<-outer(g,g,"!="); off<-row(dmat)!=col(dmat); mean(dmat[gm & off])/mean(dmat[!gm & off])}
  obs<-pseudoF(stage,D); perms<-expand.grid(rep(list(0:1),9)); pv<-numeric(nrow(perms))
  for(k in seq_len(nrow(perms))){g<-stage; for(i in 1:9) if(perms[k,i]==1) g[pair==unique(pair)[i]]<-1-g[pair==unique(pair)[i]]; pv[k]<-pseudoF(g,D)}
  beta<-data.frame(marker=marker,test="exact_paired_stage",pseudoF=obs,p=(1+sum(pv>=obs))/(1+length(pv)),permutations=length(pv))
  loo<-do.call(rbind,lapply(unique(sm$greenhouse_id),function(drop){keep<-sm$greenhouse_id!=drop; data.frame(marker=marker,dropped_greenhouse=drop,pseudoF=pseudoF(stage[keep],D[keep,keep,drop=FALSE]))}))
  genus<-if("Genus"%in%names(tx)) tx$Genus else tx$genus; genus[is.na(genus)|genus==""]<-"Unclassified"
  gx<-rowsum(ct,genus,reorder=FALSE); gx<-t(gx[,sm$sample_id,drop=FALSE])
  cs<-c(.1,.5,1); target<-cand$genus[cand$marker==marker]; rob<-do.call(rbind,lapply(cs,function(pc){clr<-log(gx+pc)-rowMeans(log(gx+pc)); do.call(rbind,lapply(target,function(g){if(!g%in%colnames(clr)) return(data.frame(marker=marker,genus=g,pseudocount=pc,mean_delta=NA,median_delta=NA,pairs_positive=NA)); delta<-sapply(unique(pair),function(p) clr[sm$sample_id==paste0("A",p),g]-clr[sm$sample_id==paste0("M",p),g]); data.frame(marker=marker,genus=g,pseudocount=pc,mean_delta=mean(delta),median_delta=median(delta),pairs_positive=sum(delta>0))}))}))
  write.table(beta,file.path(out,paste0(marker,"_exact_paired_beta.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(loo,file.path(out,paste0(marker,"_leave_one_greenhouse_out.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(rob,file.path(out,paste0(marker,"_candidate_clr_sensitivity.tsv")),sep="\t",row.names=FALSE,quote=FALSE)
  invisible(list(rel=rel,sm=sm,D=D))
}
r16<-analyze_marker("16S",file.path(base,"extension_analysis_20260822/input/16S_clean_asv_count_table.tsv"),file.path(base,"extension_analysis_20260822/input/16S_clean_taxonomy.tsv"))
rits<-analyze_marker("ITS",file.path(base,"results/dada2_ITS_v1/asv_count_table.tsv"),file.path(base,"results/its_unite_v2/its_asv_taxonomy.tsv"))
num<-setdiff(names(soil),c("sample_id","phase","lab_sample_id")); pairs<-paste0(rep(2:4,each=3),"-",1:3)
sd<-do.call(rbind,lapply(pairs,function(p){a<-soil[soil$sample_id==paste0("A",p),num,drop=FALSE]; m<-soil[soil$sample_id==paste0("M",p),num,drop=FALSE]; data.frame(pair_id=p,greenhouse_id=paste0("GH",substr(p,1,1)),a-m)})); names(sd)[-(1:2)]<-num
pc<-prcomp(sd[,num],scale.=TRUE); write.table(data.frame(sd[,1:2],PC1=pc$x[,1],PC2=pc$x[,2]),file.path(out,"soil_delta_PCA_scores.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(data.frame(variable=rownames(pc$rotation),PC1=pc$rotation[,1],PC2=pc$rotation[,2]),file.path(out,"soil_delta_PCA_loadings.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
pvfiles<-Sys.glob(file.path(base,"public_validation/results/final_cross_project_evidence_v3/*.tsv")); if(length(pvfiles)) file.copy(pvfiles,out,overwrite=TRUE)
writeLines(c("Morchella extension statistics completed","Design: 9 paired positions nested in 3 greenhouses","Beta P values use all 2^9 within-pair label swaps","Leave-one-greenhouse-out results are sensitivity estimates, not independent tests","CLR candidate results span pseudocounts 0.1, 0.5 and 1"),file.path(out,"STATISTICAL_EXTENSION_REPORT.txt"))
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))

