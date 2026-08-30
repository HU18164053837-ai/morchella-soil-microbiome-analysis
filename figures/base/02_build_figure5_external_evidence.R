options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args))args[1] else "C:/morchella_analysis"
out<-file.path(root,"manuscript_v1/figures"); src<-file.path(out,"source_data")
rd<-function(x)read.delim(x,check.names=FALSE,quote="",comment.char="")
wt<-function(x,n)write.table(x,file.path(src,n),sep="\t",quote=FALSE,row.names=FALSE,na="")
evdir<-file.path(root,"public_validation/results/final_cross_project_evidence_v3")
ef<-rd(file.path(evdir,"final_external_effect_sizes.tsv")); mx<-rd(file.path(evdir,"final_cross_project_evidence_matrix.tsv"))
gorder<-rev(c("Morchella","Mortierella","Alternaria","Botryotrichum","Terrimonas","Gemmata","Chitinophaga","Pseudarthrobacter"))
ef$genus<-factor(ef$genus,gorder)
ef$contrast_short<-ifelse(grepl("bare",ef$tested_contrast),"System 1: bare soil",ifelse(grepl("conidial",ef$tested_contrast),"System 1: conidial",""))
ef$contrast_short[grepl("primordium",ef$tested_contrast)]<-"System 1: primordium"
ef$contrast_short[ef$tested_contrast=="M1_vs_CK1"]<-"System 2: M1/CK1"
ef$contrast_short[ef$tested_contrast=="M2_vs_CK2"]<-"System 2: M2/CK2"
ef$contrast_short<-factor(ef$contrast_short,c("System 1: bare soil","System 1: conidial","System 1: primordium","System 2: M1/CK1","System 2: M2/CK2"))
status_col<-c(direction_supported="#2A9D8F",direction_opposed="#D45D4C",detected_no_clear_effect="#B7B7B7",not_detected="#ECECEC",not_testable="#FFFFFF")
grade_col<-c(externally_replicated="#117864",partially_supported="#5AB4AC",context_dependent="#E6A23C",externally_opposed="#C44E52",insufficient_external_evidence="#BDBDBD")
theme_pub<-theme_classic(base_size=7.3,base_family="Arial")+theme(axis.line=element_line(linewidth=.35),axis.ticks=element_line(linewidth=.35),legend.title=element_blank(),legend.position="top",plot.title=element_text(face="bold",size=8),plot.subtitle=element_text(size=6.8,colour="#555555"),plot.margin=margin(3,4,3,3))
panel_tag<-function(p,t)p+labs(tag=t)+theme(plot.tag=element_text(face="bold",size=9),plot.tag.position=c(.01,.99))

status_lab<-c(direction_supported="Supported",direction_opposed="Opposed",detected_no_clear_effect="No clear effect",not_detected="Not detected",not_testable="Not testable")
p5a<-ggplot(ef,aes(contrast_short,genus,fill=direction_class))+geom_tile(colour="white",linewidth=.7)+scale_fill_manual(values=status_col,labels=status_lab,drop=FALSE)+labs(x=NULL,y=NULL,title="Direction across external contrasts")+theme_pub+theme(axis.text.x=element_text(angle=35,hjust=1),legend.position="bottom",legend.text=element_text(size=5.6),legend.key.width=unit(6,"pt"),legend.spacing.x=unit(2,"pt"))+guides(fill=guide_legend(nrow=2,byrow=TRUE))

ef$system<-ifelse(grepl("System 1",ef$contrast_short),"External system 1","External system 2")
ef$plot_log2FC<-pmax(-10,pmin(10,ef$log2FC))
p5b<-ggplot(ef,aes(plot_log2FC,genus,colour=direction_class,shape=system))+geom_vline(xintercept=0,linetype=2,colour="#888888",linewidth=.35)+geom_point(size=2.0,alpha=.9,position=position_jitter(height=.11,width=0))+scale_colour_manual(values=status_col,drop=FALSE,guide="none")+scale_shape_manual(values=c("External system 1"=16,"External system 2"=17),guide="none")+coord_cartesian(xlim=c(-10,10))+labs(x="External log2 fold change",y=NULL,title="Effect-size distribution",subtitle="Circle: system 1; triangle: system 2")+theme_pub+theme(legend.position="none")

mx$genus<-factor(mx$genus,gorder)
p5c<-ggplot(mx,aes(1,genus,fill=final_evidence_grade))+geom_tile(colour="white",linewidth=.7)+geom_text(aes(label=gsub("_","\n",final_evidence_grade)),size=2.05,lineheight=.85)+scale_fill_manual(values=grade_col)+scale_x_continuous(NULL,breaks=NULL)+labs(y=NULL,title="Evidence grade")+theme_pub+theme(legend.position="none",axis.line=element_blank(),axis.ticks=element_blank())
wt(ef,"Figure5ab_external_effects.tsv");wt(mx,"Figure5c_final_evidence_grades.tsv")

draw<-function(){grid.newpage();pushViewport(viewport(layout=grid.layout(1,3,widths=unit(c(1.35,1.05,.8),"null"))));print(panel_tag(p5a,"a"),vp=viewport(layout.pos.col=1));print(panel_tag(p5b,"b"),vp=viewport(layout.pos.col=2));print(panel_tag(p5c,"c"),vp=viewport(layout.pos.col=3))}
stem<-file.path(out,"Figure5_external_validation_evidence");w<-183/25.4;h<-92/25.4
svglite::svglite(paste0(stem,".svg"),width=w,height=h);draw();dev.off()
grDevices::cairo_pdf(paste0(stem,".pdf"),width=w,height=h,family="Arial");draw();dev.off()
tmp<-file.path(Sys.getenv("TEMP"),"morchella_tiff");dir.create(tmp,recursive=TRUE,showWarnings=FALSE)
ragg::agg_tiff(file.path(tmp,"Figure5_external_validation_evidence.tiff"),width=w,height=h,units="in",res=600,compression="lzw");draw();dev.off()
ragg::agg_png(file.path(tmp,"Figure5_external_validation_evidence.png"),width=w,height=h,units="in",res=300);draw();dev.off()
writeLines(capture.output(sessionInfo()),file.path(out,"Figure5_sessionInfo.txt"))
message("FIGURE5_COMPLETE")

