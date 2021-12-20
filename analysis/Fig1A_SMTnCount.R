source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')
source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')
setwd('/data/LyuLin/Spatial/Workdirectory/SMT_figures/')

mmSI<-CreateSeuratFromSMT('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT01/ST2008PTPT01/outs',taxa.level="species")
mmSI<-getPartition(mmSI,partition.num=3)
mmSI1<-subsetPartition(mmSI,1)
mmSI2<-subsetPartition(mmSI,2)
mmSI3<-subsetPartition(mmSI,3)

SI1<-SMTPlot(mmSI1,"nCount_Spatial",max.cutoff = 'q99',colors = c("grey","orange","red"),pt.size.factor=2,plot.margin=0,legend.margin=0)
SI2<-SMTPlot(mmSI2,"nCount_Spatial",max.cutoff = 'q99',colors = c("grey","orange","red"),pt.size.factor=2,plot.margin=0,legend.margin=0)
SI3<-SMTPlot(mmSI3,"nCount_Spatial",max.cutoff = 'q99',colors = c("grey","orange","red"),pt.size.factor=2,plot.margin=0,legend.margin=0)

mmLI<-CreateSeuratFromSMT('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT03/ST2008PTPT03/outs',taxa.level="genus")
mmLI<-getPartition(mmLI,partition.num=3,MinPts=5)
mmLI1<-subsetPartition(mmLI,1)
mmLI2<-subsetPartition(mmLI,2)
mmLI3<-subsetPartition(mmLI,3)
LI1<-SMTPlot(mmLI1,"nCount_Spatial",max.cutoff = 'q99',colors = c("grey","orange","red"),pt.size.factor=2,plot.margin=0,legend.margin=0)
LI2<-SMTPlot(mmLI2,"nCount_Spatial",max.cutoff = 'q99',colors = c("grey","orange","red"),pt.size.factor=2,plot.margin=0,legend.margin=0)
LI3<-SMTPlot(mmLI3,"nCount_Spatial",max.cutoff = 'q99',colors = c("grey","orange","red"),pt.size.factor=2,plot.margin=0,legend.margin=0)

hsP53A<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR1X/outs',taxa.level="genus")
hsP53A<-getPartition(hsP53A,partition.num=2)
P531<-SMTPlot(hsP53A,"nCount_Spatial",max.cutoff = 'q99',plot.margin=0,legend.margin=0,colors=c("grey","orange","red"),pt.size.factor = 1)

hsP53B<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR2X/outs',taxa.level="genus")
hsP53B<-getPartition(hsP53B,partition.num=2)
P532<-SMTPlot(hsP53B,"nCount_Spatial",max.cutoff = 'q99',plot.margin=0,legend.margin=0,colors=c("grey","orange","red"),pt.size.factor = 1)

hsP52C<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR3X/outs',taxa.level="genus")
hsP52C<-getPartition(hsP52C,partition.num=2)
P521<-SMTPlot(hsP52C,"nCount_Spatial",max.cutoff = 'q99',plot.margin=0,legend.margin=0,colors=c("grey","orange","red"),pt.size.factor = 1)

hsP52D<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR4X/outs',taxa.level="genus")
hsP52D<-getPartition(hsP52D,partition.num=2)
P522<-SMTPlot(hsP52D,"nCount_Spatial",max.cutoff = 'q99',plot.margin=0,legend.margin=0,colors=c("grey","orange","red"),pt.size.factor = 1)

SMTPlot(hsP53A,"nCount_Spatial",pt.size.factor=0,plot.margin=0)+NoLegend()+
  SMTPlotMod(hsP53A,"nCount_Spatial",title = NULL,pt.size = 0.7)
ggsave('Fig1Ahs1.pdf',width=6.76,height=3.15)
hsP53Asubset<-readRDS('hsP53Asubset.rds')
SMTPlot(hsP53Asubset,"nCount_Spatial",pt.size.factor=0,plot.margin=0)+NoLegend()+
  SMTPlotMod(hsP53Asubset,"nCount_Spatial",title = NULL,alpha = c(0.5,1))
ggsave('Fig1Ahs1subset.pdf',width=6.76,height=3.15)
SMTPlot(hsP52C,"nCount_Spatial",pt.size.factor=0,plot.margin=0)+NoLegend()+
  SMTPlotMod(hsP52C,"nCount_Spatial",pt.size = 0.7,alpha = c(1,1),title = NULL)
ggsave('Fig1Ahs2.pdf',width=6.76,height=3.15)
SMTPlot(mmSI2,"nCount_Spatial",pt.size.factor=0,plot.margin=0)+NoLegend()+
  SMTPlotMod(mmSI2,"nCount_Spatial",title = NULL,pt.size = 2)
ggsave('Fig1Amm1.pdf',width = 6.76,height = 2.66)
SMTPlot(mmLI1,"nCount_Spatial",pt.size.factor=0,plot.margin=0)+NoLegend()+
  SMTPlotMod(mmLI1,"nCount_Spatial",title = NULL,pt.size = 2)
ggsave('Fig1Amm2.pdf',width = 6.76,height = 3.08)

