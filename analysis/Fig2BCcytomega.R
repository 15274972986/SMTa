source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')
source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')
setwd('/data/LyuLin/Spatial/Workdirectory/SMT_figures/')

library(ggrepel)
library(DESeq2)

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

hsP53<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR2X/outs',taxa.level="genus")
hsP53host<-Load10X_Spatial('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR2X/outs',filter.matrix=F,filename="raw_feature_bc_matrix.h5")

SMTPlot(hsP53,"genus-Cytomegalovirus",legend.title="Cytomegalovirus",pt.size.factor=0,colors=c("blue","orange","red"))+NoLegend()
SMTPlotMod(hsP53,"genus-Cytomegalovirus",title = "Cytomegalovirus")
ggsave('Fig2B1.pdf',width=5.6,height=5.6)
targetspots=filterMicrobeSpot(hsP53,feature="genus-Cytomegalovirus",min.count=8)
viewSpots(hsP53,spots=targetspots)

hsP53host<-AddModuleScore(hsP53host,features=list("Fibroblast"=c("LUM","DCN","VIM","PDGFRA","COL1A2")))
SMTPlotMod(hsP53host,"Cluster1",title = "fibroblast score")
ggsave('Fig2B2.pdf',width=5.6,height=5.6)
fibroblastRichSpots=Cells(subset(hsP53host,Cluster1>4))
viewSpots(hsP53host,spots=fibroblastRichSpots,color = "blue")
viewSpots(hsP53,spots=targetspots)

virusContainingSpots=intersect(targetspots,fibroblastRichSpots)
virusAbsentSpots=setdiff(fibroblastRichSpots,targetspots)

df<-hsP53host@assays$Spatial@counts[,c(virusContainingSpots,virusAbsentSpots)] %>% as.data.frame()
grouplist<-c(rep("virus_enriched",length(virusContainingSpots)),rep("virus_free",length(virusAbsentSpots)))
MetaData<-data.frame(row.names=colnames(df), group_list=grouplist)
df<-df[(df %>% rowSums())!=0,]
dds<-DESeqDataSetFromMatrix(countData=df,colData=MetaData,design=~group_list)
dds2<-DESeq(dds)
res<-results(dds2, contrast=c("group_list","virus_enriched","virus_free"))
res<-as.data.frame(res)
res$gene<-rownames(res)
genelist=(arrange(res,by=pvalue) %>% dplyr::filter(.,pvalue<=0.01))$gene
#GeneEnrich(genelist,"hs")
#ggsave("Fig2C.pdf",width=10.5,height=4.62)
data=arrange(res,by=pvalue) %>% dplyr::filter(.,pvalue<=0.0005)
res[,"color"]<-ifelse(res$pvalue<0.01&abs(res$log2FoldChange)>log2(1.5),"red","grey")
ggplot(res)+geom_point(aes(x=log2FoldChange,y=-log10(pvalue),color=color))+scale_color_manual(values = c("grey","red"))+
  geom_hline(yintercept=2,color="red",linetype="longdash")+xlim(c(-3,3))+ylim(c(0,6))+geom_vline(xintercept=log2(1.5),color="red",linetype="longdash")+
  geom_vline(xintercept=-log2(1.5),color="red",linetype="longdash")+
  geom_text_repel(data=res %>% dplyr::filter(.,gene %in% c("HLA-E","TNIP1","IL32")),aes(x=log2FoldChange,y=-log10(pvalue),label=gene))+
  geom_point(data=res %>% dplyr::filter(.,gene %in% c("HLA-E","TNIP1","IL32")),aes(x=log2FoldChange,y=-log10(pvalue)),color="blue")+
  theme_minimal(base_size = 15)+NoLegend()
ggsave('Fig2C.pdf',width=4.22,height=5.64)
