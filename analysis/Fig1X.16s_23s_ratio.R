source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')
source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')
setwd('/data/LyuLin/Spatial/Workdirectory/SMT_figures/')

data0<-read.delim('16s23s.tsv',header = T)
colnames(data0)[3:7]<-c("16s.count","23s.count","16s.ratio","23s.ratio","16s.23s.ratio")
data0<-data0[,1:4]
data0$others<-data0$Total.bacterial.read.count-data0$`16s.count`-data0$`23s.count`
rownames(data0)<-data0$sample
data0<-data0[,3:5] %>% t %>% as.data.frame()
data<-TransformRawDf(data0)
data<-data %>% group_by(feat1)
total<-summarise(data,total=sum(values))
data<-left_join(data,total,by="feat1")
data$values<-data$values/data$total
data$feat1<-c(rep("human",12),rep("mouse",12))
data<-data %>% group_by(feat1,feat2)
data<-summarise(data,ratio=mean(values))
colnames(data)[3]<-"count"

plot_ratio<-function(data,sample){
  data=dplyr::filter(data,feat1==sample)
  colnames(data)=c("sample","origin","ratio")
  data$origin=data$origin %>% gsub(".count","",.,fixed = T)
  data$origin=paste0(data$origin," (",round(data$ratio*100/sum(data$ratio),2),"%)")
  ggplot(data)+geom_bar(aes(x=sample,y=ratio,fill=origin),stat="identity",position ='stack')+
    scale_fill_simpsons()+coord_polar(theta='y')+theme_void()+
    theme(legend.margin=margin(r=0.2,unit="cm"),text=element_text(size=15),
          title=element_text(size=10),legend.title=element_text(size=12))+
    ggtitle(sample)
}

plots<-list()
samples<-data$feat1 %>% unique()
i<-1
for(sample in samples){
  plots[[i]]=plot_ratio(data,sample)
  i=i+1
}
ggarrange(plotlist=plots,nrow=1,ncol=2)
ggsave('Fig1D_16s23sRatio.pdf',width = 6.51,height=1.69)
