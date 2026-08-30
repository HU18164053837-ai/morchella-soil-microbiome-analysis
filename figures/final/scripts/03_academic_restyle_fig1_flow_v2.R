options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(ggplot2); library(grid); library(svglite); library(ragg); library(scales)})

# Figure contract
# Core conclusion: morel cultivation is associated with coordinated bacterial–fungal
# community restructuring, whereas individual genera show heterogeneous external transferability.
# Archetypes: Fig. 1 schematic-led composite; Figs. 2–6 quantitative/asymmetric grids.
# Integrity: all source rows are retained; no inferential statistic is recomputed here.

args <- commandArgs(trailingOnly=TRUE)
manu <- if(length(args)) gsub("\\\\","/",args[1]) else stop("manuscript_v1 path required")
old <- file.path(manu,"figures","source_data")
old2 <- file.path(manu,"figures_v2","source_data")
out <- if(length(args) >= 2) gsub("\\\\","/",args[2]) else file.path(manu,"figures_academic_v2_20260823")
src <- file.path(out,"source_data")
dir.create(out,recursive=TRUE,showWarnings=FALSE); dir.create(src,recursive=TRUE,showWarnings=FALSE)

rd <- function(path) read.delim(path,check.names=FALSE,quote="",comment.char="")
copy_src <- function(path) {file.copy(path,file.path(src,basename(path)),overwrite=TRUE); rd(path)}

COL <- c(
  mist="#FCE8E6", coral="#FFC6BC", rose="#F8B9B8",
  cloud="#D6DFEF", ice="#A5CDE2", lake="#5FA3CB", sky="#668FCA",
  lilac="#C9CEFE", sea="#66A8CD", deep="#006BAD",
  ink="#263238", grey="#69757D", pale="#EDF2F6", white="#FFFFFF"
)
stage_col <- c(M="#F8B9B8", A="#006BAD")
gh_col <- c(GH2="#A5CDE2",GH3="#668FCA",GH4="#006BAD")
domain_col <- c(`16S`="#006BAD", ITS="#7D82D6")
direction_col <- c(Decrease="#F8B9B8",Increase="#006BAD")

theme_pub <- function(base=6.8) theme_classic(base_size=base,base_family="Arial") + theme(
  axis.line=element_line(linewidth=.32,colour=COL["ink"]),axis.ticks=element_line(linewidth=.32,colour=COL["ink"]),
  axis.text=element_text(colour=COL["ink"]),axis.title=element_text(colour=COL["ink"]),
  legend.title=element_blank(),legend.text=element_text(size=5.8),legend.key.height=unit(7,"pt"),
  strip.background=element_blank(),strip.text=element_text(face="bold",colour=COL["ink"]),
  plot.title=element_text(face="bold",size=7.3,colour=COL["ink"],hjust=0),
  plot.subtitle=element_text(size=5.8,colour=COL["grey"],margin=margin(b=3)),
  panel.grid=element_blank(),plot.margin=margin(3,4,3,4))
theme_set(theme_pub())

tag <- function(p,x) p + labs(tag=x) + theme(plot.tag=element_text(face="bold",size=8.0,colour=COL["ink"]),plot.tag.position=c(.005,.995))

draw_layout <- function(items,nrow,ncol,widths=rep(1,ncol),heights=rep(1,nrow)) {
  grid.newpage(); pushViewport(viewport(layout=grid.layout(nrow,ncol,widths=unit(widths,"null"),heights=unit(heights,"null"))))
  for(it in items) print(it$plot,vp=viewport(layout.pos.row=it$row,layout.pos.col=it$col))
}
save_multi <- function(name,items,nrow,ncol,widths=rep(1,ncol),heights=rep(1,nrow),width_mm=169,height_mm=120) {
  w<-width_mm/25.4; h<-height_mm/25.4; stem<-file.path(out,name)
  raster_dir <- file.path(Sys.getenv("TEMP"),"ase_v3_rasters"); dir.create(raster_dir,recursive=TRUE,showWarnings=FALSE)
  svglite::svglite(paste0(stem,".svg"),width=w,height=h); draw_layout(items,nrow,ncol,widths,heights); dev.off()
  grDevices::cairo_pdf(paste0(stem,".pdf"),width=w,height=h,family="Arial"); draw_layout(items,nrow,ncol,widths,heights); dev.off()
  ragg::agg_tiff(file.path(raster_dir,paste0(name,".tiff")),width=w,height=h,units="in",res=600,compression="lzw"); draw_layout(items,nrow,ncol,widths,heights); dev.off()
  ragg::agg_png(file.path(raster_dir,paste0(name,".png")),width=w,height=h,units="in",res=300); draw_layout(items,nrow,ncol,widths,heights); dev.off()
}

# ---------- Figure 1: schematic-led design and QC ----------
meta <- copy_src(file.path(old,"Figure1a_sample_design.tsv"))
trim <- copy_src(file.path(old,"Figure1c_primer_retention.tsv"))
final <- copy_src(file.path(old,"Figure1d_final_data.tsv"))
flow <- data.frame(
  x=1:6,
  phase=factor(c("Sampling","Sampling","Amplicon profiling","Amplicon profiling","Community analysis","Community analysis"),
               levels=c("Sampling","Amplicon profiling","Community analysis")),
  label=c("Paired soil\nsampling","24 soils\n48 libraries","Full-length\n16S + ITS","Primer-aware\nquality control","DADA2 ASVs\nSILVA / UNITE","Cross-domain +\nexternal evidence")
)
flow$xmin <- flow$x - 0.39; flow$xmax <- flow$x + 0.39
phase_cols <- c("Sampling"=unname(COL["rose"]),"Amplicon profiling"=unname(COL["lake"]),"Community analysis"=unname(COL["deep"]))
arrow_df <- data.frame(x=1.40:5.40, xend=1.60:5.60)
phase_df <- data.frame(x=c(1.5,3.5,5.5), label=c("Sampling design","Amplicon profiling","Community analysis"))
p1a <- ggplot() +
  geom_segment(data=arrow_df,aes(x=x,xend=xend,y=0.53,yend=0.53),
               arrow=arrow(length=unit(1.25,"mm"),type="closed"),linewidth=.38,colour=COL["grey"])+
  geom_rect(data=flow,aes(xmin=xmin,xmax=xmax,ymin=.27,ymax=.78),
            fill="white",colour=COL["cloud"],linewidth=.48)+
  geom_rect(data=flow,aes(xmin=xmin,xmax=xmax,ymin=.73,ymax=.78,fill=phase),colour=NA)+
  geom_text(data=flow,aes(x=xmin+.07,y=.67,label=sprintf("%02d",x),colour=phase),
            hjust=0,size=2.05,fontface="bold")+
  geom_text(data=flow,aes(x=x,y=.48,label=label),size=2.05,lineheight=.92,colour=COL["ink"])+
  geom_segment(data=data.frame(x=c(.61,2.61,4.61),xend=c(2.39,4.39,6.39)),
               aes(x=x,xend=xend,y=.93,yend=.93),linewidth=.42,colour=COL["cloud"])+
  geom_text(data=phase_df,aes(x=x,y=.98,label=label),size=2.05,fontface="bold",colour=COL["grey"])+
  scale_fill_manual(values=phase_cols,guide="none")+scale_colour_manual(values=phase_cols,guide="none")+
  coord_cartesian(xlim=c(.48,6.52),ylim=c(.20,1.04),clip="off")+
  theme_void(base_family="Arial")+theme(plot.margin=margin(4,10,1,10))

paired <- meta[meta$cohort=="cultivation_paired",]
paired$stage <- ifelse(paired$condition=="before_cultivation","M","A")
paired$x <- as.numeric(factor(paired$pair_id,levels=unique(paired$pair_id)))
p1b <- ggplot(paired,aes(x,stage,group=pair_id))+
  geom_segment(aes(xend=x,y="M",yend="A"),linewidth=1.7,colour=COL["cloud"],lineend="round")+
  geom_point(aes(fill=stage,colour=greenhouse_id),shape=21,size=3.0,stroke=.8)+
  scale_fill_manual(values=stage_col)+scale_colour_manual(values=gh_col)+
  scale_x_continuous(breaks=1:9,labels=paired$pair_id[paired$stage=="M"],expand=expansion(mult=.06))+
  labs(x="Matched spatial position",y=NULL,title="Sampling pairs",subtitle=NULL)+
  theme_pub()+theme(legend.position="none",axis.line.y=element_blank(),axis.ticks.y=element_blank())

p1c <- ggplot(trim,aes(retention,marker,colour=marker))+
  geom_segment(aes(x=98.8,xend=retention,yend=marker),linewidth=5,colour=COL["pale"],lineend="round")+
  geom_point(size=4)+geom_text(aes(label=sprintf("%.2f%%",retention)),hjust=-.35,size=2.7,fontface="bold",colour=COL["ink"])+
  scale_colour_manual(values=domain_col,guide="none")+scale_x_continuous(limits=c(98.8,100.2),breaks=c(99,99.5,100))+
  labs(x="Reads retained after primer trimming",y=NULL,title="Primer retention",subtitle=NULL)+
  theme_pub()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())

cards <- rbind(transform(final,metric="ASVs",value=ASVs),transform(final,metric="Reads",value=reads))
cards$metric <- factor(cards$metric,c("ASVs","Reads"))
p1d <- ggplot(cards,aes(metric,value,fill=marker))+
  geom_col(width=.58,position=position_dodge(.68),colour=COL["white"],linewidth=.3)+
  geom_text(aes(label=comma(value)),position=position_dodge(.68),vjust=-.45,size=2.45,fontface="bold",colour=COL["ink"])+
  scale_fill_manual(values=domain_col)+scale_y_log10(labels=label_number(scale_cut=cut_si("")),expand=expansion(mult=c(.05,.18)))+
  labs(x=NULL,y="Final count (log scale)",title="Final feature recovery",subtitle=NULL)+
  theme_pub()+theme(legend.position="top",legend.direction="horizontal",legend.background=element_blank())
save_multi("Figure1_design_and_data_quality_v3",list(list(plot=tag(p1a,"a"),row=1,col=1:3),list(plot=tag(p1b,"b"),row=2,col=1),list(plot=tag(p1c,"c"),row=2,col=2),list(plot=tag(p1d,"d"),row=2,col=3)),2,3,widths=c(1.35,1,1),heights=c(.8,1),height_mm=120)

# ---------- Figure 2: beta restructuring as hero evidence ----------
b16 <- copy_src(file.path(old,"Figure2a_16S_PCoA.tsv")); fits <- copy_src(file.path(old,"Figure2b_ITS_PCoA.tsv"))
a16 <- copy_src(file.path(old,"Figure2c_16S_Shannon.tsv")); aits <- copy_src(file.path(old,"Figure2d_ITS_Shannon_deltas.tsv"))
pcoa <- function(d,title,sub){
  ggplot(d,aes(axis1,axis2,group=pair_id))+
    geom_path(colour=COL["cloud"],linewidth=.55,arrow=arrow(length=unit(1.5,"mm"),type="closed"))+
    geom_point(aes(fill=stage,colour=greenhouse_id),shape=21,size=3.2,stroke=.75)+
    scale_fill_manual(values=stage_col)+scale_colour_manual(values=gh_col)+
    labs(x=sprintf("PCoA1 (%.1f%%)",d$axis1_percent[1]),y=sprintf("PCoA2 (%.1f%%)",d$axis2_percent[1]),title=title,subtitle=sub)+
    theme_pub(7.5)+theme(aspect.ratio=1,legend.position="top",legend.box="horizontal",legend.spacing.x=unit(2,"pt"))
}
p2a <- pcoa(b16,"16S community structure","Bray-Curtis R2 = 0.152; exploratory P = 0.0039")
p2b <- pcoa(fits,"ITS community structure","Bray-Curtis R2 = 0.149; exploratory P = 0.0117")

a16$stage<-factor(a16$stage,c("M","A"))
p2c <- ggplot(a16,aes(stage,shannon,group=pair_id))+
  geom_line(aes(colour=greenhouse_id),linewidth=.65,alpha=.72)+geom_point(aes(fill=stage,colour=greenhouse_id),shape=21,size=2.65,stroke=.65)+
  scale_fill_manual(values=stage_col)+scale_colour_manual(values=gh_col)+
  labs(x=NULL,y="Shannon diversity",title="16S Shannon diversity",subtitle="Mean A - M = +0.206; 6/9 pairs increased")+
  theme_pub()+theme(legend.position="none")
q <- aits[aits$subset=="UNITE_kingdom_Fungi",]
ql <- rbind(data.frame(pair_id=q$pair_id,greenhouse_id=q$greenhouse_id,stage="M",shannon=q$shannon_before),data.frame(pair_id=q$pair_id,greenhouse_id=q$greenhouse_id,stage="A",shannon=q$shannon_after)); ql$stage<-factor(ql$stage,c("M","A"))
p2d <- ggplot(ql,aes(stage,shannon,group=pair_id))+
  geom_line(aes(colour=greenhouse_id),linewidth=.65,alpha=.72)+geom_point(aes(fill=stage,colour=greenhouse_id),shape=21,size=2.65,stroke=.65)+
  scale_fill_manual(values=stage_col)+scale_colour_manual(values=gh_col)+
  labs(x=NULL,y="Shannon diversity",title="ITS Shannon diversity",subtitle="Mean A - M = +0.942; 6/9 pairs increased")+
  theme_pub()+theme(legend.position="none")
save_multi("Figure2_paired_community_restructuring_v3",list(list(plot=tag(p2a,"a"),row=1,col=1),list(plot=tag(p2b,"b"),row=1,col=2),list(plot=tag(p2c,"c"),row=2,col=1),list(plot=tag(p2d,"d"),row=2,col=2)),2,2,heights=c(1.35,.8),height_mm=145)

# ---------- Figure 3: sample-resolved composition heatmaps + directional genus effects ----------
bp <- copy_src(file.path(old2,"Figure3a_16S_phylum_composition.tsv")); fp <- copy_src(file.path(old2,"Figure3b_ITS_phylum_composition.tsv"))
bg <- copy_src(file.path(old2,"Figure3c_16S_abundance_weighted_genus_shifts.tsv")); fg <- copy_src(file.path(old2,"Figure3d_ITS_abundance_weighted_genus_shifts.tsv"))
sample_level <- function(ids){
  num <- as.numeric(gsub("[^0-9]","",sub("-.*","",ids))); pos<-as.numeric(sub(".*-","",ids)); st<-substr(ids,1,1)
  order(num,pos,factor(st,c("M","A")))
}
heat <- function(d,taxon,title){
  d$sample_id<-factor(d$sample_id,levels=unique(d$sample_id[sample_level(d$sample_id)])); d[[taxon]]<-factor(d[[taxon]],levels=rev(unique(d[[taxon]])))
  ggplot(d,aes(sample_id,.data[[taxon]],fill=relative_abundance))+
    geom_tile(colour="white",linewidth=.22)+
    scale_fill_gradientn(colours=c(COL["mist"],COL["cloud"],COL["ice"],COL["lake"],COL["deep"]),values=rescale(c(0,.01,.05,.2,1)),labels=percent_format(accuracy=1),name="Relative\nabundance")+
    labs(x="Matched samples ordered M to A within position",y=NULL,title=title,subtitle="Samples ordered by matched position")+
    theme_pub()+theme(axis.text.x=element_text(angle=55,hjust=1,size=5.4),legend.position="right",legend.title=element_text(size=6,face="bold"),axis.line=element_blank(),axis.ticks=element_blank())
}
p3a<-heat(bp,"taxon","16S phylum composition")
p3b<-heat(fp,"taxon_plot","ITS phylum composition")
effect <- function(d,title){
  d$direction<-ifelse(d$effect>0,"Increase","Decrease"); d$genus<-factor(d$genus,levels=d$genus[order(d$effect)])
  ggplot(d,aes(effect,genus))+
    geom_vline(xintercept=0,colour=COL["grey"],linetype=2,linewidth=.35)+
    geom_segment(aes(x=0,xend=effect,yend=genus,colour=direction),linewidth=1.15,lineend="round")+
    geom_point(aes(size=mean_abundance,fill=direction),shape=21,colour="white",stroke=.45)+
    scale_colour_manual(values=direction_col,guide="none")+scale_fill_manual(values=direction_col)+scale_size_continuous(range=c(2,5.8),labels=percent_format(accuracy=.1))+
    labs(x="log2 mean relative-abundance ratio (A/M)",y=NULL,title=title,subtitle="Point size: larger stage mean")+
    theme_pub()+theme(legend.position="bottom",legend.box="horizontal")
}
p3c<-effect(bg,"16S genus-level shifts")
p3d<-effect(fg,"ITS genus-level shifts")
save_multi("Figure3_community_composition_and_genus_shifts_v3",list(list(plot=tag(p3a,"a"),row=1,col=1),list(plot=tag(p3b,"b"),row=1,col=2),list(plot=tag(p3c,"c"),row=2,col=1),list(plot=tag(p3d,"d"),row=2,col=2)),2,2,heights=c(1.18,.82),height_mm=150)

# ---------- Figure 4: cross-domain coordination ----------
bd<-copy_src(file.path(old,"Figure3a_cross_marker_Bray_distances.tsv")); ad<-copy_src(file.path(old,"Figure3bc_cross_marker_alpha_deltas.tsv"))
p4a<-ggplot(bd,aes(distance_16S,distance_ITS))+
  geom_point(shape=21,fill=COL["ice"],colour=COL["deep"],size=1.55,stroke=.35,alpha=.62)+
  geom_smooth(method="lm",se=TRUE,linewidth=1,colour=COL["deep"],fill=COL["cloud"],alpha=.55)+
  annotate("label",x=Inf,y=-Inf,label="Mantel r = 0.590\nP = 0.0001",hjust=1.12,vjust=-.45,size=2.65,label.size=0,fill=COL["mist"],colour=COL["ink"])+
  labs(x="Bacterial Bray-Curtis distance",y="Fungal Bray-Curtis distance",title="Cross-domain beta diversity",subtitle="All pairwise distances among the 18 matched-stage samples")+
  theme_pub(7.8)
quad <- function(x,y,xlab,ylab,title,sub){
  ggplot(ad,aes(.data[[x]],.data[[y]]))+
    annotate("rect",xmin=-Inf,xmax=0,ymin=-Inf,ymax=0,fill=COL["mist"],alpha=.4)+annotate("rect",xmin=0,xmax=Inf,ymin=0,ymax=Inf,fill=COL["cloud"],alpha=.45)+
    geom_hline(yintercept=0,colour=COL["grey"],linewidth=.35)+geom_vline(xintercept=0,colour=COL["grey"],linewidth=.35)+
    geom_point(aes(fill=greenhouse_id),shape=21,size=3.2,colour="white",stroke=.5)+scale_fill_manual(values=gh_col)+
    labs(x=xlab,y=ylab,title=title,subtitle=sub)+theme_pub()+theme(legend.position="bottom")
}
p4b<-quad("shannon_delta_16S","shannon_delta_ITS","Bacterial Delta Shannon","Fungal Delta Shannon","Shannon diversity changes","Spearman rho = -0.167; P = 0.668")
p4c<-quad("observed_delta_16S","observed_delta_ITS","Bacterial Delta Observed ASVs","Fungal Delta Observed ASVs","Observed-ASV changes","Spearman rho = -0.159; P = 0.683")
save_multi("Figure4_cross_domain_coordination_v3",list(list(plot=tag(p4a,"a"),row=1,col=1),list(plot=tag(p4b,"b"),row=1,col=2),list(plot=tag(p4c,"c"),row=1,col=3)),1,3,widths=c(1.35,1,1),height_mm=78)

# ---------- Figure 5: integrated frozen-candidate evidence portrait ----------
cg<-copy_src(file.path(old,"Figure4_frozen_candidate_summary.tsv")); ra<-copy_src(file.path(old,"Figure4b_candidate_mean_relative_abundance.tsv"))
ord<-rev(c("Morchella","Mortierella","Alternaria","Botryotrichum","Terrimonas","Gemmata","Chitinophaga","Pseudarthrobacter")); cg$genus<-factor(cg$genus,ord); ra$genus<-factor(ra$genus,ord)
cg$direction<-ifelse(cg$log2FC>0,"Increase","Decrease")
p5a<-ggplot(cg,aes(log2FC,genus))+
  geom_vline(xintercept=0,linetype=2,colour=COL["grey"],linewidth=.35)+geom_segment(aes(x=0,xend=log2FC,yend=genus,colour=direction),linewidth=1.2,lineend="round")+
  geom_point(aes(fill=marker,shape=tier),size=3.1,colour="white",stroke=.5)+scale_colour_manual(values=direction_col,guide="none")+scale_fill_manual(values=domain_col)+scale_shape_manual(values=c(Primary=21,Supporting=24))+
  labs(x="Local log2 ratio (A/M)",y=NULL,title="Local effect size",subtitle="Fill: marker; triangle: supporting")+theme_pub()+theme(legend.position="none")
p5b<-ggplot(ra,aes(plot_abundance,genus,group=genus))+
  geom_line(colour=COL["cloud"],linewidth=1.2)+geom_point(aes(fill=stage),shape=21,size=3.2,colour="white",stroke=.5)+scale_fill_manual(values=stage_col)+
  scale_x_log10(limits=c(.001,40),breaks=c(.001,.01,.1,1,10),labels=function(x)paste0(x,"%"))+
  labs(x="Mean relative abundance (log scale)",y=NULL,title="Mean relative abundance",subtitle="Display floor; values unchanged")+theme_pub()+theme(legend.position="bottom")
p5c<-ggplot(cg,aes(direction_pct,genus,fill=tier))+
  geom_col(width=.52,colour="white",linewidth=.35)+geom_text(aes(label=paste0(pmax(pairs_up,pairs_down),"/9")),hjust=1.15,size=2.45,fontface="bold",colour="white")+
  scale_fill_manual(values=c(Primary="#006BAD",Supporting="#C9CEFE"))+scale_x_continuous(limits=c(0,105),breaks=c(0,50,100),labels=function(x)paste0(x,"%"))+
  labs(x="Pairs matching direction",y=NULL,title="Directional consistency",subtitle="Greenhouse means reviewed")+theme_pub()+theme(legend.position="bottom")
save_multi("Figure5_frozen_candidate_genera_v3",list(list(plot=tag(p5a,"a"),row=1,col=1),list(plot=tag(p5b,"b"),row=1,col=2),list(plot=tag(p5c,"c"),row=1,col=3)),1,3,widths=c(1.12,1.12,1),height_mm=100)

# ---------- Figure 6: external stress-test evidence ----------
ef<-copy_src(file.path(old,"Figure5ab_external_effects.tsv")); mx<-copy_src(file.path(old,"Figure5c_final_evidence_grades.tsv")); ef$genus<-factor(ef$genus,ord); mx$genus<-factor(mx$genus,ord)
ef$contrast_short<-factor(ef$contrast_short,levels=c("System 1: bare soil","System 1: conidial","System 1: primordium","System 2: M1/CK1","System 2: M2/CK2"))
status_col<-c(direction_supported="#006BAD",direction_opposed="#F8B9B8",detected_no_clear_effect="#D6DFEF",not_detected="#F3F5F6",not_testable="#FFFFFF")
status_lab<-c(direction_supported="Supports",direction_opposed="Opposes",detected_no_clear_effect="Unclear",not_detected="Not detected",not_testable="Not testable")
p6a<-ggplot(ef,aes(contrast_short,genus,fill=direction_class))+
  geom_tile(colour="white",linewidth=1)+scale_fill_manual(values=status_col,labels=status_lab,drop=FALSE)+
  labs(x=NULL,y=NULL,title="External direction concordance",subtitle="Independent public datasets")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))+theme_pub()+theme(axis.text.x=element_text(angle=32,hjust=1,size=5.7),legend.position="bottom",axis.line=element_blank(),axis.ticks=element_blank())
ef$system2<-ifelse(grepl("System 1",ef$contrast_short),"System 1","System 2")
p6b<-ggplot(ef,aes(plot_log2FC,genus))+
  geom_vline(xintercept=0,linetype=2,colour=COL["grey"],linewidth=.35)+
  geom_segment(aes(x=0,xend=plot_log2FC,yend=genus,colour=direction_class),linewidth=.65,alpha=.55)+
  geom_point(aes(fill=direction_class,shape=system2),size=2.35,colour="white",stroke=.4)+
  scale_colour_manual(values=status_col,guide="none")+scale_fill_manual(values=status_col,guide="none")+
  scale_shape_manual(values=c("System 1"=21,"System 2"=24),labels=c("System 1 (circle)","System 2 (triangle)"))+
  coord_cartesian(xlim=c(-10,10))+labs(x="External log2 ratio (clipped at +/-10)",y=NULL,title="External effect sizes",subtitle="Circles: system 1; triangles: system 2")+
  theme_pub()+theme(legend.position="none")
grade_col<-c(partially_supported="#668FCA",context_dependent="#C9CEFE",externally_opposed="#F8B9B8",insufficient_external_evidence="#D6DFEF",externally_replicated="#006BAD")
grade_lab<-c(partially_supported="Partial\nsupport",context_dependent="Context\ndependent",externally_opposed="Externally\nopposed",insufficient_external_evidence="Insufficient\nevidence",externally_replicated="Replicated")
p6c<-ggplot(mx,aes(1,genus,fill=final_evidence_grade))+
  geom_tile(colour="white",linewidth=1)+geom_text(aes(label=grade_lab[final_evidence_grade]),size=2.3,lineheight=.82,colour=COL["ink"])+
  scale_fill_manual(values=grade_col)+scale_x_continuous(NULL,breaks=NULL)+labs(y=NULL,title="Evidence grade",subtitle="Discordant evidence retained")+
  theme_pub()+theme(legend.position="none",axis.line=element_blank(),axis.ticks=element_blank())
save_multi("Figure6_external_validation_evidence_v3",list(list(plot=tag(p6a,"a"),row=1,col=1),list(plot=tag(p6b,"b"),row=1,col=2),list(plot=tag(p6c,"c"),row=1,col=3)),1,3,widths=c(1.28,1.04,.8),height_mm=112)

# ---------- Supplementary Figure S1 ----------
soil<-copy_src(file.path(old2,"SupplementaryFigureS1a_soil_paired_changes.tsv")); cor<-copy_src(file.path(old2,"SupplementaryFigureS1b_candidate_soil_correlations.tsv")); hd16<-copy_src(file.path(old2,"SupplementaryFigureS1c_16S_HD_alpha.tsv")); hdits<-copy_src(file.path(old2,"SupplementaryFigureS1d_ITS_HD_alpha.tsv"))
soil$metric<-factor(soil$metric,levels=unique(soil$metric)); cor$genus<-factor(cor$genus,ord)
pS1a<-ggplot(soil,aes(standardized_change,metric,colour=greenhouse_id))+
  geom_vline(xintercept=0,colour=COL["grey"],linewidth=.35)+geom_segment(aes(x=0,xend=standardized_change,yend=metric),linewidth=.55,alpha=.45)+geom_point(size=2.45)+
  scale_colour_manual(values=gh_col)+labs(x="Standardized paired change (A - M)",y=NULL,title="Soil properties",subtitle="Nine positions; raw units retained in source data")+theme_pub()+theme(legend.position="bottom")
pS1b<-ggplot(cor,aes(metric,genus,fill=spearman_rho))+
  geom_tile(colour="white",linewidth=.65)+geom_text(aes(label=sprintf("%.2f",spearman_rho)),size=2.15,colour=COL["ink"])+
  scale_fill_gradient2(low=unname(COL["rose"]),mid=unname(COL["mist"]),high=unname(COL["deep"]),midpoint=0,limits=c(-1,1),name="Spearman rho")+
  labs(x=NULL,y=NULL,title="Candidate-soil correlations",subtitle="All BH-adjusted q > 0.05; exploratory only")+theme_pub()+theme(axis.text.x=element_text(angle=42,hjust=1,size=5.8),axis.line=element_blank(),axis.ticks=element_blank(),legend.position="bottom")
raw_alpha<-function(d,title){
  ggplot(d,aes(condition,value,fill=condition))+
    geom_point(shape=21,size=2.7,colour="white",stroke=.45,position=position_jitter(width=.09,height=0))+
    stat_summary(fun=mean,geom="crossbar",width=.48,linewidth=.55,colour=COL["ink"])+facet_wrap(~metric,scales="free_y")+
    scale_fill_manual(values=c("Healthy point"="#A5CDE2","Diseased point"="#F8B9B8"))+
    labs(x=NULL,y=NULL,title=title,subtitle="One point per condition; three spatial subsamples")+theme_pub()+theme(axis.text.x=element_text(angle=25,hjust=1,size=5.6),legend.position="none")
}
pS1c<-raw_alpha(hd16,"Bacterial alpha diversity")
pS1d<-raw_alpha(hdits,"Fungal alpha diversity")
save_multi("SupplementaryFigureS1_soil_and_white_mold_context_v3",list(list(plot=tag(pS1a,"a"),row=1,col=1),list(plot=tag(pS1b,"b"),row=1,col=2),list(plot=tag(pS1c,"c"),row=2,col=1),list(plot=tag(pS1d,"d"),row=2,col=2)),2,2,widths=c(1,1.18),heights=c(1.12,.88),height_mm=145)

writeLines(c(capture.output(sessionInfo()),"",paste("Source rows were copied without exclusion to",src)),file.path(out,"sessionInfo.txt"))
cat("REDESIGN_COMPLETE\n")







