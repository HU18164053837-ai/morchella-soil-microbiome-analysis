options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "C:/morchella_analysis"
out <- file.path(root,"manuscript_v1/figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src <- file.path(out,"source_data")
dir.create(src,recursive=TRUE,showWarnings=FALSE)

rd <- function(x) read.delim(x,check.names=FALSE,quote="",comment.char="")
wt <- function(x,n) write.table(x,file.path(src,n),sep="\t",quote=FALSE,row.names=FALSE,na="")

pal <- c(M="#3B6FB6",A="#D66A4A",ITS="#8A6BBE",S16="#3B8C88",neutral="#6E6E6E",light="#D9D9D9",core="#315A9A",support="#C28B35")
pv <- function(x) unname(pal[x])
theme_pub <- function(base=7.4) theme_classic(base_size=base,base_family="Arial") + theme(
  axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),
  legend.title=element_blank(),legend.position="top",legend.key.width=unit(8,"pt"),
  strip.background=element_blank(),strip.text=element_text(face="bold"),
  plot.title=element_text(face="bold",size=base+.4,hjust=0),
  plot.subtitle=element_text(size=base-.4,colour="#555555"),panel.grid=element_blank(),
  plot.margin=margin(3,4,3,3))
theme_set(theme_pub())

panel_tag <- function(p,tag) p + labs(tag=tag) + theme(plot.tag=element_text(face="bold",size=9),plot.tag.position=c(.01,.99))
draw_grid <- function(plots,nrow,ncol,widths=rep(1,ncol),heights=rep(1,nrow)) {
  grid.newpage(); pushViewport(viewport(layout=grid.layout(nrow,ncol,widths=unit(widths,"null"),heights=unit(heights,"null"))))
  for(i in seq_along(plots)){r <- ((i-1)%/%ncol)+1; c <- ((i-1)%%ncol)+1; print(plots[[i]],vp=viewport(layout.pos.row=r,layout.pos.col=c))}
}
save_grid <- function(stem,plots,nrow,ncol,width_mm=183,height_mm=125,widths=rep(1,ncol),heights=rep(1,nrow)) {
  w <- width_mm/25.4; h <- height_mm/25.4
  svglite::svglite(paste0(stem,".svg"),width=w,height=h); draw_grid(plots,nrow,ncol,widths,heights); dev.off()
  grDevices::cairo_pdf(paste0(stem,".pdf"),width=w,height=h,family="Arial"); draw_grid(plots,nrow,ncol,widths,heights); dev.off()
  tiff_tmp_dir <- file.path(Sys.getenv("TEMP"),"morchella_tiff")
  dir.create(tiff_tmp_dir,recursive=TRUE,showWarnings=FALSE)
  tiff_tmp <- file.path(tiff_tmp_dir,paste0(basename(stem),".tiff"))
  ragg::agg_tiff(tiff_tmp,width=w,height=h,units="in",res=600,compression="lzw"); draw_grid(plots,nrow,ncol,widths,heights); dev.off()
  png_tmp <- file.path(tiff_tmp_dir,paste0(basename(stem),".png"))
  ragg::agg_png(png_tmp,width=w,height=h,units="in",res=300); draw_grid(plots,nrow,ncol,widths,heights); dev.off()
}

# Figure 1: design and QC
meta <- rd(file.path(root,"metadata/sample_metadata.tsv"))
meta_bio <- unique(meta[,c("biological_sample_id","cohort","condition","continuous_cropping_years","greenhouse_id","spatial_point","pair_id")])
trim <- rd(file.path(root,"qc/primer_trimming_summary.tsv"))
trim_sum <- aggregate(cbind(input_reads,trimmed_reads)~marker,trim,sum)
trim_sum$retention <- trim_sum$trimmed_reads/trim_sum$input_reads*100
qc_final <- data.frame(marker=c("16S","ITS"),ASVs=c(9477,8056),reads=c(188611,942492))
wt(meta_bio,"Figure1a_sample_design.tsv"); wt(trim_sum,"Figure1c_primer_retention.tsv"); wt(qc_final,"Figure1d_final_data.tsv")

design <- data.frame(x=c(1,2,3,4,5,6),y=1,label=c("Soil\nsampling","24 biological\nsamples","Full-length\n16S + ITS","Primer\ntrimming","DADA2 +\ntaxonomy","Community\nintegration"))
p1a <- ggplot(design,aes(x,y)) + geom_segment(aes(x=1,xend=6,y=1,yend=1),colour=pal["light"],linewidth=1.2) +
  geom_point(shape=21,size=5,fill=c("#E9EEF7","#E9EEF7","#DDEDEB","#F5E6DF","#E9EEF7","#E5DFF0"),colour="#444444",stroke=.45) +
  geom_text(aes(label=label),vjust=-1.55,size=2.0,lineheight=.88) + coord_cartesian(xlim=c(.7,6.3),ylim=c(.75,1.38),clip="off") + theme_void(base_family="Arial") + theme(plot.margin=margin(10,8,5,8))

hier <- data.frame(group=c("M/A paired","H/D case"),greenhouses=c(3,1),independent_sites=c(9,2),spatial_subsamples=c(18,6))
hier_long <- reshape(hier,varying=c("greenhouses","independent_sites","spatial_subsamples"),v.names="count",timevar="level",times=c("Greenhouses","Sampling sites","Soil samples"),direction="long")
hier_long$level <- factor(hier_long$level,c("Greenhouses","Sampling sites","Soil samples"))
p1b <- ggplot(hier_long,aes(level,count,fill=group)) + geom_col(position=position_dodge(.7),width=.62,colour="white") +
  geom_text(aes(label=count),position=position_dodge(.7),vjust=-.25,size=2.4) + scale_fill_manual(values=c("M/A paired"=pv("S16"),"H/D case"=pv("ITS"))) +
  scale_y_continuous(expand=expansion(mult=c(0,.16))) + labs(x=NULL,y="Count",subtitle="Sampling hierarchy (subsamples are not independent units)")

p1c <- ggplot(trim_sum,aes(marker,retention,fill=marker)) + geom_segment(aes(xend=marker,y=99,yend=retention),linewidth=2.2,colour="#D9D9D9") + geom_point(shape=21,size=3.2,colour="white",stroke=.4) + geom_text(aes(label=sprintf("%.2f%%",retention)),vjust=-.9,size=2.6) +
  scale_fill_manual(values=c("16S"=pv("S16"),"ITS"=pv("ITS"))) + scale_y_continuous(limits=c(98.8,100.05),breaks=c(99,99.5,100),expand=c(0,0)) + labs(x=NULL,y="Reads retained (%)",subtitle="Primer orientation and trimming") + theme(legend.position="none")

qcl <- rbind(data.frame(marker=qc_final$marker,metric="ASVs",value=qc_final$ASVs,label=format(qc_final$ASVs,big.mark=",")),data.frame(marker=qc_final$marker,metric="Reads",value=qc_final$reads,label=format(qc_final$reads,big.mark=",")))
p1d <- ggplot(qcl,aes(marker,value,fill=marker)) + geom_col(width=.58) + geom_text(aes(label=label),vjust=-.35,size=2.35) + facet_wrap(~metric,scales="free_y") +
  scale_fill_manual(values=c("16S"=pv("S16"),"ITS"=pv("ITS"))) + scale_y_continuous(labels=scales::label_comma(),expand=expansion(mult=c(0,.18))) + labs(x=NULL,y="Final non-chimeric data",subtitle="All 24 samples retained") + theme(legend.position="none")
save_grid(file.path(out,"Figure1_design_and_data_quality"),list(panel_tag(p1a,"a"),panel_tag(p1b,"b"),panel_tag(p1c,"c"),panel_tag(p1d,"d")),2,2,height_mm=126,widths=c(1.35,1),heights=c(.8,1))

# Figure 2: paired community shifts
s16p <- rd(file.path(root,"downstream/16S_v1_silva_diversity/beta/pcoa_coordinates.tsv")); s16p <- s16p[s16p$metric=="Bray_Curtis" & s16p$cohort=="cultivation_paired",]
itsp <- rd(file.path(root,"downstream/ITS_v2_taxonomy/diversity/fungal_bray_pcoa.tsv")); itsp <- itsp[itsp$cohort=="cultivation_paired",]
stage <- function(x) ifelse(x=="before_cultivation","M","A")
s16p$stage <- factor(stage(s16p$condition),c("M","A")); itsp$stage <- factor(stage(itsp$condition),c("M","A"))
wt(s16p,"Figure2a_16S_PCoA.tsv"); wt(itsp,"Figure2b_ITS_PCoA.tsv")
pcoa_plot <- function(d,title,sub,col) ggplot(d,aes(axis1,axis2,group=pair_id)) + geom_path(colour="#B7B7B7",linewidth=.35,arrow=arrow(length=unit(1.2,"mm"),type="closed")) + geom_point(aes(fill=stage),shape=21,size=2.4,colour="white",stroke=.35) + scale_fill_manual(values=c(M=pv("M"),A=pv("A"))) + labs(x=sprintf("PCoA1 (%.1f%%)",d$axis1_percent[1]),y=sprintf("PCoA2 (%.1f%%)",d$axis2_percent[1]),title=title,subtitle=sub) + coord_equal() + theme(legend.position="top")
p2a <- pcoa_plot(s16p,"Bacterial community","Bray-Curtis R2=0.152; exploratory P=0.0039")
p2b <- pcoa_plot(itsp,"Fungal community","Bray-Curtis R2=0.149; exploratory P=0.0117")

a16 <- rd(file.path(root,"downstream/16S_v1_silva_diversity/alpha/alpha_diversity_depth1500.tsv")); a16 <- a16[a16$cohort=="cultivation_paired",]; a16$stage <- factor(stage(a16$condition),c("M","A"))
aits_delta <- rd(file.path(root,"downstream/ITS_v3_design_aware/alpha/MA_paired_changes_all_vs_fungal.tsv"))
wt(a16[,c("sample_id","pair_id","greenhouse_id","stage","shannon")],"Figure2c_16S_Shannon.tsv"); wt(aits_delta,"Figure2d_ITS_Shannon_deltas.tsv")
p2c <- ggplot(a16,aes(stage,shannon,group=pair_id)) + geom_line(colour="#AFAFAF",linewidth=.4) + geom_point(aes(fill=stage),shape=21,size=2.1,colour="white",stroke=.3) + scale_fill_manual(values=c(M=pv("M"),A=pv("A"))) + labs(x=NULL,y="Shannon diversity",title="Bacterial alpha diversity",subtitle="Mean paired change = +0.206; 6/9 increased") + theme(legend.position="none")
if(all(c("subset","pair_id","shannon_change") %in% names(aits_delta))){q <- aits_delta[aits_delta$subset=="UNITE_kingdom_Fungi",]; p2d <- ggplot(q,aes(reorder(pair_id,shannon_change),shannon_change,fill=shannon_change>0)) + geom_col(width=.65) + geom_hline(yintercept=0,linewidth=.35) + coord_flip() + scale_fill_manual(values=c(`TRUE`=pv("A"),`FALSE`=pv("M"))) + labs(x="Paired sampling point",y="A - M Shannon",title="Fungal alpha diversity",subtitle="Mean paired change = +0.942; 6/9 increased") + theme(legend.position="none")} else {p2d <- ggplot(data.frame(x=1,y=1),aes(x,y))+annotate("text",1,1,label="ITS Shannon source table schema requires review")+theme_void()}
save_grid(file.path(out,"Figure2_paired_community_restructuring"),list(panel_tag(p2a,"a"),panel_tag(p2b,"b"),panel_tag(p2c,"c"),panel_tag(p2d,"d")),2,2,height_mm=145)

# Figure 3: beta coordination and alpha divergence
bray <- function(mat){mat<-as.matrix(mat); mat<-sweep(mat,2,colSums(mat),"/"); n<-ncol(mat); z<-matrix(0,n,n,dimnames=list(colnames(mat),colnames(mat))); for(i in 1:(n-1)) for(j in (i+1):n){z[i,j]<-z[j,i]<-sum(abs(mat[,i]-mat[,j]))/2}; z}
c16 <- rd(file.path(root,"downstream/16S_v1_silva_diversity/objects/16S_clean_asv_count_table.tsv")); rownames(c16)<-c16$ASV_ID; c16$ASV_ID<-NULL
cits <- rd(file.path(root,"downstream/ITS_v2_taxonomy/tables/fungal_asv_count_table.tsv")); rownames(cits)<-cits$ASV_ID; cits$ASV_ID<-NULL
ma <- intersect(colnames(c16),colnames(cits)); ma <- ma[grepl("^[MA]",ma)]
d16 <- bray(c16[,ma,drop=FALSE]); dits <- bray(cits[,ma,drop=FALSE]); ix <- lower.tri(d16)
bd <- data.frame(distance_16S=d16[ix],distance_ITS=dits[ix]); wt(bd,"Figure3a_cross_marker_Bray_distances.tsv")
p3a <- ggplot(bd,aes(distance_16S,distance_ITS)) + geom_point(size=1.2,alpha=.55,colour=pv("neutral")) + geom_smooth(method="lm",se=TRUE,linewidth=.7,colour=pv("core"),fill="#CBD8EA") + labs(x="16S Bray-Curtis distance",y="ITS Bray-Curtis distance",title="Cross-domain beta structure",subtitle="Mantel r=0.590; P=0.0001 (M/A samples)")
ad <- rd(file.path(root,"downstream/joint_16S_ITS_v1/tables/MA_cross_marker_alpha_changes.tsv")); wt(ad,"Figure3bc_cross_marker_alpha_deltas.tsv")
p3b <- ggplot(ad,aes(shannon_delta_16S,shannon_delta_ITS)) + geom_hline(yintercept=0,colour=pv("light"),linewidth=.35) + geom_vline(xintercept=0,colour=pv("light"),linewidth=.35) + geom_point(aes(fill=greenhouse_id),shape=21,size=2.6,colour="white",stroke=.35) + geom_smooth(method="lm",se=FALSE,linewidth=.6,colour=pv("neutral")) + scale_fill_manual(values=c(GH2="#5B8DB8",GH3="#5A9F73",GH4="#B87979")) + labs(x="16S Delta Shannon (A - M)",y="ITS Delta Shannon (A - M)",title="Shannon responses",subtitle="Spearman rho=-0.167; P=0.668")
p3c <- ggplot(ad,aes(observed_delta_16S,observed_delta_ITS)) + geom_hline(yintercept=0,colour=pv("light"),linewidth=.35) + geom_vline(xintercept=0,colour=pv("light"),linewidth=.35) + geom_point(aes(fill=greenhouse_id),shape=21,size=2.6,colour="white",stroke=.35) + geom_smooth(method="lm",se=FALSE,linewidth=.6,colour=pv("neutral")) + scale_fill_manual(values=c(GH2="#5B8DB8",GH3="#5A9F73",GH4="#B87979")) + labs(x="16S Delta Observed ASVs (A - M)",y="ITS Delta Observed ASVs (A - M)",title="Richness responses",subtitle="Spearman rho=-0.159; P=0.683")
save_grid(file.path(out,"Figure3_cross_domain_coordination"),list(panel_tag(p3a,"a"),panel_tag(p3b,"b"),panel_tag(p3c,"c")),1,3,height_mm=68)

# Figure 4: pre-frozen genera
cg <- rd(file.path(root,"core_genera_audit_20260822/eight_genera_audit_summary.tsv")); cg$genus <- factor(cg$genus,rev(c("Morchella","Mortierella","Alternaria","Botryotrichum","Terrimonas","Gemmata","Chitinophaga","Pseudarthrobacter")))
cg$log2FC <- log2((cg$mean_A_RA+1e-5)/(cg$mean_M_RA+1e-5)); cg$direction_pct <- pmax(cg$pairs_up,cg$pairs_down)/9*100; cg$tier <- ifelse(cg$freeze_status=="PRIMARY_CORE","Primary","Supporting")
wt(cg,"Figure4_frozen_candidate_summary.tsv")
p4a <- ggplot(cg,aes(log2FC,genus,colour=marker,shape=tier)) + geom_vline(xintercept=0,linewidth=.4,colour="#777777") + geom_segment(aes(x=0,xend=log2FC,yend=genus),linewidth=.7) + geom_point(size=2.8,fill="white") + scale_colour_manual(values=c(ITS=pv("ITS"),`16S`=pv("S16"))) + scale_shape_manual(values=c(Primary=16,Supporting=17)) + labs(x="log2 fold change (A/M; pseudocount 1e-5)",y=NULL,title="Signed compositional shift")
rl <- rbind(data.frame(genus=cg$genus,stage="M",abundance=cg$mean_M_RA*100,marker=cg$marker),data.frame(genus=cg$genus,stage="A",abundance=cg$mean_A_RA*100,marker=cg$marker)); rl$stage<-factor(rl$stage,c("M","A")); rl$plot_abundance<-pmax(rl$abundance,.001); wt(rl,"Figure4b_candidate_mean_relative_abundance.tsv")
p4b <- ggplot(rl,aes(plot_abundance,genus,group=genus)) + geom_line(colour="#B8B8B8",linewidth=.6) + geom_point(aes(fill=stage),shape=21,size=2.6,colour="white",stroke=.35) + scale_fill_manual(values=c(M=pv("M"),A=pv("A"))) + scale_x_log10(limits=c(.001,40),breaks=c(.001,.01,.1,1,10),labels=function(x) paste0(x,"%")) + labs(x="Mean relative abundance (log scale)",y=NULL,title="Mean abundance in M and A")
p4c <- ggplot(cg,aes(direction_pct,genus,fill=tier)) + geom_col(width=.58) + geom_text(aes(label=paste0(pmax(pairs_up,pairs_down),"/9")),hjust=1.15,colour="white",size=2.35,fontface="bold") + scale_fill_manual(values=c(Primary=pv("core"),Supporting=pv("support"))) + scale_x_continuous(limits=c(0,105),breaks=c(0,50,100),labels=function(x)paste0(x,"%")) + labs(x="Paired-direction consistency",y=NULL,title="Within-study consistency")
save_grid(file.path(out,"Figure4_frozen_candidate_genera"),list(panel_tag(p4a,"a"),panel_tag(p4b,"b"),panel_tag(p4c,"c")),1,3,height_mm=82,widths=c(1.1,1.05,.85))

writeLines(c(capture.output(sessionInfo()),"","Input root:",root),file.path(out,"figure_build_sessionInfo.txt"))
message("FIGURES_1_TO_4_COMPLETE")

