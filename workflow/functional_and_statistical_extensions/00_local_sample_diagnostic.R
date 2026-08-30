m<-read.delim("C:/morchella_analysis/metadata/sample_metadata.tsv",check.names=FALSE)
if(!"sample_id"%in%names(m)) m$sample_id<-m$biological_sample_id
for(mark in c("16S","ITS")){
 f<-if(mark=="16S") "C:/morchella_analysis/downstream/16S_v1_silva_diversity/objects/16S_clean_asv_count_table.tsv" else "C:/morchella_analysis/hpc_results/ITS_dada2_v1_corrected_v1/results/asv_count_table.tsv"
 h<-names(read.delim(f,nrows=1,check.names=FALSE)); s<-m[m$marker==mark & m$cohort=="cultivation_paired",]
 cat(mark,"meta",nrow(s),"overlap",sum(s$sample_id%in%h),"\n"); print(setdiff(s$sample_id,h))
}

