suppressPackageStartupMessages({library(ggplot2);library(grid);library(svglite);library(ragg)})

# Figure contract
# Claim: cultivation visibly restructures dominant bacterial and fungal composition,
# while the largest genus-level shifts differ between the two microbial domains.
# Archetype: quantitative grid; all observations retained in composition panels.

root <- "C:/morchella_analysis"
out <- file.path(root,"manuscript_v1","figures_v2")
src <- file.path(out,"source_data")
dir.create(out,recursive=TRUE,showWarnings=FALSE);dir.create(src,recursive=TRUE,showWarnings=FALSE)
tmp <- file.path(Sys.getenv("TEMP"),"morchella_tiff");dir.create(tmp,showWarnings=FALSE)

rtsv <- function(x) read.delim(x,check.names=FALSE,stringsAsFactors=FALSE)
wtsv <- function(x,n) write.table(x,file.path(src,n),sep="\t",row.names=FALSE,quote=FALSE)
stage_label <- function(x) ifelse(grepl("before|M_",x),"Before (M)","After (A)")
sample_order <- function(x) {
  z <- sub("^[MA]", "", x); parts <- strsplit(z,"-",fixed=TRUE)
  gh <- as.numeric(vapply(parts,function(v) v[1],"")); pos <- as.numeric(vapply(parts,function(v) v[2],""))
  order(gh,pos,substr(x,1,1)=="A")
}

pal <- c("#35618C","#4E91A8","#63B6A5","#8BC17C","#C7C66B","#E4B85C","#DD8B57","#C76563","#956B8E","#7C8795","#D9D9D9")
theme_pub <- theme_classic(base_size=7.2,base_family="Arial")+
  theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),
        legend.title=element_blank(),legend.text=element_text(size=5.8),
        plot.title=element_text(face="bold",size=8),plot.subtitle=element_text(size=6.5,colour="#555555"),
        plot.margin=margin(4,5,3,4))
tag <- function(p,x) p+labs(tag=x)+theme(plot.tag=element_text(face="bold",size=9),plot.tag.position=c(.01,.99))

# 16S phylum composition: paired M/A samples only, top taxa table already contains Other.
b <- rtsv(file.path(root,"downstream","16S_v1_silva_diversity","taxonomy","phylum_top10_long.tsv"))
b <- b[grepl("^[MA][0-9]+-[0-9]+$",b$sample_id),]
b$stage <- factor(stage_label(b$condition_short),levels=c("Before (M)","After (A)"))
b$sample_id <- factor(b$sample_id,levels=unique(b$sample_id[sample_order(b$sample_id)]))
b$taxon <- factor(b$taxon,levels=rev(unique(b$taxon[order(ave(b$relative_abundance,b$taxon,FUN=mean),decreasing=TRUE)])))

# ITS phylum composition: cultivation-paired samples; retain dominant phyla and aggregate the rest.
f <- rtsv(file.path(root,"downstream","ITS_v2_taxonomy","tables","composition_phylum.tsv"))
f <- f[f$cohort=="cultivation_paired",]
fm <- aggregate(relative_abundance~taxon,f,mean)
keepf <- head(fm$taxon[order(fm$relative_abundance,decreasing=TRUE)],8)
f$taxon_plot <- ifelse(f$taxon %in% keepf,f$taxon,"Other")
f <- aggregate(relative_abundance~sample_id+condition+taxon_plot,f,sum)
f$stage <- factor(stage_label(f$condition),levels=c("Before (M)","After (A)"))
f$sample_id <- factor(f$sample_id,levels=unique(f$sample_id[sample_order(f$sample_id)]))
f$taxon_plot <- factor(f$taxon_plot,levels=rev(c(keepf,"Other")))

comp_plot <- function(d,taxon_col,title){
  ggplot(d,aes(sample_id,relative_abundance,fill=.data[[taxon_col]]))+
    geom_col(width=.88,colour=NA)+facet_grid(~stage,scales="free_x",space="free_x")+
    scale_y_continuous("Relative abundance",labels=scales::percent_format(accuracy=1),expand=c(0,0))+
    scale_fill_manual(values=setNames(rep(pal,length.out=length(levels(d[[taxon_col]]))),levels(d[[taxon_col]])))+
    labs(x=NULL,title=title)+theme_pub+
    theme(axis.text.x=element_text(angle=60,hjust=1,size=5.7),strip.background=element_blank(),strip.text=element_text(face="bold",size=6.5),legend.position="right",legend.key.height=unit(7,"pt"))
}
p3a <- comp_plot(b,"taxon","Bacterial phylum composition")
p3b <- comp_plot(f,"taxon_plot","Fungal phylum composition")

# Genus-level overview: a pre-specified abundance-weighted ranking among taxa
# observed in at least 9 M/A samples. This is descriptive, not an FDR claim.
bs <- rtsv(file.path(root,"downstream","16S_v1_silva_diversity","taxonomy","MA_genus_screening.tsv"))
bs <- bs[bs$prevalence_MA>=9 & bs$genus!="Unclassified_at_genus",]
bs$max_mean <- pmax(bs$mean_relative_M,bs$mean_relative_A)
bs$score <- abs(bs$log2_mean_relative_ratio_A_vs_M)*bs$max_mean
bs <- head(bs[order(bs$score,decreasing=TRUE),],8)
bs$domain <- "Bacteria";bs$effect <- bs$log2_mean_relative_ratio_A_vs_M;bs$mean_abundance <- bs$max_mean

fs <- rtsv(file.path(root,"downstream","ITS_v4_genus_screening","tables","MA_genus_screening.tsv"))
fs <- fs[fs$prevalence_MA>=9,]
fs$max_mean <- pmax(fs$mean_relative_before,fs$mean_relative_after)
fs$score <- abs(fs$log2_mean_relative_ratio)*fs$max_mean
fs <- head(fs[order(fs$score,decreasing=TRUE),],8)
fs$domain <- "Fungi";fs$effect <- fs$log2_mean_relative_ratio;fs$mean_abundance <- fs$max_mean

effect_plot <- function(d,title){
  d$genus <- factor(d$genus,levels=d$genus[order(d$effect)])
  ggplot(d,aes(effect,genus,size=mean_abundance,colour=effect>0))+
    geom_vline(xintercept=0,linetype=2,colour="#888888",linewidth=.35)+geom_point(alpha=.9)+
    scale_colour_manual(values=c(`TRUE`="#D17B49",`FALSE`="#3E7EA8"),labels=c(`TRUE`="Higher after",`FALSE`="Lower after"))+
    scale_size_continuous(range=c(1.8,5.2),labels=scales::percent_format(accuracy=.1))+
    labs(x="log2 mean abundance ratio (A/M)",y=NULL,title=title,subtitle="Point size: larger group mean relative abundance")+
    theme_pub+theme(legend.position="bottom",legend.box="vertical",legend.key.width=unit(8,"pt"))+
    guides(colour=guide_legend(order=1),size=guide_legend(order=2,title="Mean abundance"))
}
p3c <- effect_plot(bs,"Abundance-weighted bacterial genus shifts")
p3d <- effect_plot(fs,"Abundance-weighted fungal genus shifts")

wtsv(b,"Figure3a_16S_phylum_composition.tsv");wtsv(f,"Figure3b_ITS_phylum_composition.tsv")
wtsv(bs,"Figure3c_16S_abundance_weighted_genus_shifts.tsv");wtsv(fs,"Figure3d_ITS_abundance_weighted_genus_shifts.tsv")

draw <- function(){
  grid.newpage();pushViewport(viewport(layout=grid.layout(2,2,widths=unit(c(1,1),"null"),heights=unit(c(1.12,.88),"null"))))
  print(tag(p3a,"a"),vp=viewport(layout.pos.row=1,layout.pos.col=1));print(tag(p3b,"b"),vp=viewport(layout.pos.row=1,layout.pos.col=2))
  print(tag(p3c,"c"),vp=viewport(layout.pos.row=2,layout.pos.col=1));print(tag(p3d,"d"),vp=viewport(layout.pos.row=2,layout.pos.col=2))
}
stem <- file.path(out,"Figure3_community_composition_and_genus_shifts");width_mm<-183;height_mm<-145;w<-width_mm/25.4;h<-height_mm/25.4
svglite::svglite(paste0(stem,".svg"),width=w,height=h);draw();dev.off()
grDevices::cairo_pdf(paste0(stem,".pdf"),width=w,height=h,family="Arial");draw();dev.off()
ragg::agg_tiff(file.path(tmp,"Figure3_community_composition_and_genus_shifts.tiff"),width=w,height=h,units="in",res=600);draw();dev.off()
ragg::agg_png(file.path(tmp,"Figure3_community_composition_and_genus_shifts.png"),width=w,height=h,units="in",res=300);draw();dev.off()
writeLines(capture.output(sessionInfo()),file.path(out,"Figure3_sessionInfo.txt"))
cat("FIGURE3_COMPLETE\n")

