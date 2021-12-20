source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')
source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')
setwd('/data/LyuLin/Spatial/Workdirectory/SMT_figures/')

getMeta<-function(srt,sample.name){
  meta=dplyr::filter(srt@meta.data,nCount_Spatial>0)
  meta$sample=sample.name
  meta=meta %>% as_tibble()
  return(meta)
}

mmSI<-CreateSeuratFromSMT('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT01/ST2008PTPT01/outs',taxa.level="genus")
mmSI<-getMeta(mmSI,"mmSI")
mmLI<-CreateSeuratFromSMT('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT03/ST2008PTPT03/outs',taxa.level="genus")
mmLI<-getMeta(mmLI,"mmLI")
hsP531<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR1X/outs',taxa.level="genus")
hsP531<-getMeta(hsP531,"hsP53-1")
hsP532<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR2X/outs',taxa.level="genus")
hsP532<-getMeta(hsP532,"hsP53-2")
hsP521<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR3X/outs',taxa.level="genus")
hsP521<-getMeta(hsP521,"hsP52-1")
hsP522<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR4X/outs',taxa.level="genus")
hsP522<-getMeta(hsP522,"hsP52-2")
data<-Reduce(rbind,list(mmSI,mmLI,hsP531,hsP532,hsP521,hsP522))

data$nCount_Rank<-case_when(data$nCount_Spatial<=10~"UMI<=10",
                               data$nCount_Spatial>10 & data$nCount_Spatial<=100~"10<UMI<=100",
                               data$nCount_Spatial>100~"UMI>100")

data<-group_by()
ggplot(data)+geom_density(aes(x=log10(nCount_Spatial),color=sample),show.legend=FALSE)+
  theme_minimal()+scale_color_d3()+xlab("log10(UMI count)")+
  stat_density(aes(x=log10(nCount_Spatial),colour=sample),geom="line",position="identity")

ggsave("Fig1E.pdf",width=8.31,height=2.59)
  

