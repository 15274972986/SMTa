source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')
source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')
setwd('/data/LyuLin/Spatial/Workdirectory/SMT_figures/')


#msrt: SMT seurat object
#hsrt: host seurat object
#assay: "default" or "Smooth"
#min.microbe: min microbe taxa number a spot should contain, otherwise it would be removed
#min.gene: min host gene number a spot should contain, otherwise it would be removed
#min.microbe.freq: min times a microbe taxa should occur in sample, otherwise it would be removed
#min.gene.freq: min times a gene should occur in sample, otherwise it would be removed
#taxa.level: taxa level that would like to find correlations
#corthreshold: R^2 value that were going to be preserved
#demo.corplus.threshold & demo.corminus.threshold: number of correlation pairs to be present in dotplot
#to get more spot, decrease
SMTCCA<-function(msrt,hsrt,normalize.by.spot=F,log.trans=T,corthreshold=0.4,demo.corplus.threshold=100,demo.corminus.threshold=100){
  mmtx0=msrt@assays$Smooth@counts %>% as.data.frame()
  if(normalize.by.spot){
    mmtx=mmtx0 %>% t %>% as.array() %>% proportions(margin=1)
  }else{
    mmtx=mmtx0 %>% t
  }
  if(log.trans){
    mmtx=mmtx %>% +1 %>% log10()
  }
  hsrt=findPeriTissue(hsrt,layerstart=0,layerend=1)
  layer1spot=hsrt@assays$Diffusion$peri.layers[[3]]
  hsrt=subset(hsrt,cells=c(layer1spot,filterTissueSpots(hsrt)))
  hmtx0=hsrt@assays$Smooth@counts %>% as.data.frame()
  if(normalize.by.spot){
    hmtx=hmtx0 %>% t %>% as.array() %>% proportions(margin=1)
  }else{
    hmtx=hmtx0 %>% t
  }
  if(log.trans){
    hmtx=hmtx %>% +1 %>% log10()
  }
  sharedspot=intersect(rownames(mmtx),rownames(hmtx))
  mmtx=mmtx[sharedspot,]
  hmtx=hmtx[sharedspot,]
  message("Performing CCA ...")
  cca=matcor(mmtx,hmtx)
  raw=cca$XYcor %>% as.data.frame()
  raw=raw[colnames(mmtx),colnames(hmtx)]
  raw[is.na(raw)]=0
  strong=TransformRawDf(raw,corthreshold)
  strong$feat1=strong$feat1
  strong$feat2=strong$feat2
  strong=strong[order(strong$values,decreasing=T),]
  plotdata=rbind(strong %>% head(.,demo.corplus.threshold),strong %>% tail(.,demo.corminus.threshold))
  plotdata$correlation=ifelse(plotdata$values>0,"#DC143C","#0000FF")
  plotdata$feat2=plotdata$feat2 %>% gsub("[A-z]+-","",.)
  p=ggplot(plotdata)+geom_point(aes(feat2,feat1,size=values %>% abs(),color=correlation),color=plotdata$correlation)+
    theme_cowplot()+
    theme(text = element_text(size=12),axis.text.y=element_text(size = 12),
          plot.title = element_text(hjust = 0.5),axis.text.x=element_text(angle=45, hjust=1,face = "italic",size=12))+
    ylab("Host genes")+xlab("microbes")+scale_size_continuous(name="correlation")
  res=list()
  res[["demo"]]<-p
  res[["microbe.mtx"]]=mmtx0[colnames(mmtx),sharedspot]
  res[["host.mtx"]]=hmtx0[colnames(hmtx),sharedspot]
  res[["cca"]]=raw
  res[["strong.correlation"]]=strong
  res[["microbe.richness"]]=(res$microbe.mtx!=0) %>% table %>% proportions()
  res[["host.richness"]]=(res$host.mtx!=0) %>% table %>% proportions()
  return(res)
}

plot.pair<-function(out,index){
  feat1=out$strong.correlation$feat1[index]
  feat2=out$strong.correlation$feat2[index]
  SMTPlotMod(mmSI1,feat2,colors = c("red","red"))+SMTPlotMod(mmSIhost2,feat1,colors = c("cyan","blue"),alpha = c(0,10))
}

plot.cor<-function(out,plus,minus){
  out$cca=dplyr::select(out$cca,!starts_with("Gm"))
  out$cca=dplyr::select(out$cca,!ends_with("Rik"))
  data=out$cca %>% TransformRawDf(.,min(abs(minus),plus))
  data=data %>% dplyr::filter(.,values>=plus|values<=minus)
  data$feat2=gsub("genus-","",data$feat2)
  data$color=ifelse(data$values>0,"#DC143C","#0000FF")
  ggplot(data)+geom_point(aes(x=feat2,y=feat1,size=abs(values)),color=data$color)+
    theme_cowplot()+
    theme(text = element_text(size=12),axis.text.y=element_text(size = 12),
          plot.title = element_text(hjust = 0.5),axis.text.x=element_text(angle=45, hjust=1,face = "italic",size=12))+
    ylab("Host genes")+xlab("microbes")+scale_size_continuous(name="correlation")
}

GeneEnrich<-function(gene,species){
  if(species=="hs"){
    db="org.Hs.eg.db"
  }else if(species=="mm"){
    db="org.Mm.eg.db"
  }else{
    stop("Unsupported species")
  }
  genelist=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb=db)
  res=enrichGO(gene=genelist$ENTREZID, OrgDb=get(db), keyType = "ENTREZID", ont = c("ALL"),
      pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2,
      minGSSize = 10, maxGSSize = 500, readable = FALSE, pool = FALSE)
  dotplot(res,font.size=8)
}

mmSI<-CreateSeuratFromSMT('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT01/ST2008PTPT01/outs',taxa.level="genus")
mmSIhost<-Load10X_Spatial('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT01/ST2008PTPT01/outs',filter.matrix=F,filename="raw_feature_bc_matrix.h5")
#mmSIhost.fil<-Load10X_Spatial('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT01/ST2008PTPT01/outs',filter.matrix=T,filename="filtered_feature_bc_matrix.h5")
mmSI<-denoise(mmSI,n.permutation=1000,strict.mode = F)
mmSI<-SmoothSpots(mmSI,mode="smt")
#mmSIhost<-SCTransform(mmSIhost,assay = "Spatial")
mmSIhost<-subset(mmSIhost,features=rownames(mmSIhost)[rowSums(mmSIhost)>2])
hvg<-(FindVariableFeatures(mmSIhost,nfeatures=10000))@assays$Spatial@var.features
mmSIhost<-subset(mmSIhost,features=hvg)
mmSIhost<-SmoothSpots(mmSIhost)
SMTPlot(mmSIhost,"nCount_Spatial",colors=c("blue","orange","red"),
        pt.size.factor=1.6,legend.title = "nCount_Spatial",legend.text.size = 15)+
  SMTPlot(mmSIhost,"nCount_Smooth",colors=c("blue","orange","red"),
          pt.size.factor=1.6,legend.title = "nCount_Smooth",legend.text.size = 15)
ggsave("FigS2A.pdf",width=9.43,height=5.67)
out<-SMTCCA(mmSI,mmSIhost,corthreshold=0.3,log.trans = F,normalize.by.spot = F)
plot.cor(out,0.3,0.3)
ggsave('Fig2Amod.pdf',width = 10,height=12)
DefaultAssay(mmSI)<-"Phyloseq_level"
DefaultAssay(mmSIhost)<-"Smooth"
SMTPlotMod(mmSI,"genus-Porphyromonas",colors = c("lightgrey","red"))+SMTPlotMod(mmSIhost,"Saa1",colors = c("lightgrey","blue"),alpha = c(0,1))
ggsave('Fig2C_Porphyromonas_vs_Saa1.pdf',width = 9.48,height=6.02)

#GeneEnrich((out$strong.correlation %>% dplyr::filter(.,feat2=="genus-Lactobacillus" & values<=0))$feat1,"mm")
#ggsave('FigS2A.pdf',width=8.92,height=4.28)

mmLI<-CreateSeuratFromSMT('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT03/ST2008PTPT03/outs',taxa.level="genus")
mmLIhost<-Load10X_Spatial('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT03/ST2008PTPT03/outs',filter.matrix=F,filename="raw_feature_bc_matrix.h5")
mmLI<-denoise(mmLI,n.permutation = 1000,strict.mode = F)
mmLI<-SmoothSpots(mmLI,mode="smt")
mmLIhost<-subset(mmLIhost,features=rownames(mmLIhost)[rowSums(mmLIhost)>2])
hvg2<-(FindVariableFeatures(mmLIhost,nfeatures=10000))@assays$Spatial@var.features
mmLIhost<-subset(mmLIhost,features=hvg2)
mmLIhost<-SmoothSpots(mmLIhost)
out2<-SMTCCA(mmLI,mmLIhost,corthreshold=0.3,log.trans = F,normalize.by.spot = F)
plot.cor(out2,0.4,0.4)+coord_flip()
ggsave("FigS2BMod.pdf",width=12.8,height=2.62)
DefaultAssay(mmLI)<-"Smooth"
DefaultAssay(mmLIhost)<-"Smooth"
SMTPlotMod(mmLI,"genus-Helicobacter",colors = c("lightgrey","red"),alpha=c(0,1))+
  SMTPlotMod(mmLIhost,"Saa1",colors = c("lightgrey","blue"),alpha = c(0,1))
ggsave('Fig2D_Helicobacter_vs_Saa1.pdf',width = 9.48,height=6.02)
SMTPlotMod(mmLI,"genus-Helicobacter",colors = c("lightgrey","red"),alpha=c(0,1))+
  SMTPlotMod(mmLIhost,"Dmbt1",colors = c("lightgrey","blue"),alpha = c(0,1))
ggsave('Fig2D_Helicobacter_vs_Dmbt1.pdf',width = 9.48,height=6.02)

save.image('cca2.data')

mmSI<-getPartition(mmSI,partition.num=3)
mmSI1<-subsetPartition(mmSI,1)
mmSI2<-subsetPartition(mmSI,2)
mmSI3<-subsetPartition(mmSI,3)
mmSIhost<-getPartition(mmSIhost,partition.num=3)
mmSIhost1<-subsetPartition(mmSIhost,1)
mmSIhost2<-subsetPartition(mmSIhost,2)
mmSIhost3<-subsetPartition(mmSIhost,3)
out3<-SMTCCA(mmSI1,mmSIhost,corthreshold=0.3,log.trans = F,normalize.by.spot = F)
plot.cor(out3,0.45,0.45)
ggsave('Fig2Amod2.pdf',width = 10,height=12)
SMTPlotMod(mmSI1,"genus-Staphylococcus",colors = c("red","red"))+SMTPlotMod(mmSIhost2,"Defa17",colors = c("cyan","blue"),alpha = c(0,10))
