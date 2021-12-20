source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')
source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')

library('org.Mm.eg.db')
library('AnnotationHub')

setwd('/data/LyuLin/Spatial/Workdirectory/BulkRNA-2020-11-12/blastraw/')

coplot<-function(x,y,df){
  feat1<-x
  feat2<-y
  ggscatter(df, x = feat1, y = feat2, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = as.character(feat1), ylab = as.character(feat2),color = "#0095DE")
}

getBulkTaxaData<-function(level){
  allsample=Reduce(rbind,list(
#                   CallTaxaTableFromTSV('TS2010CLZM01.merged.tsv',level,"TS2010CLZM01"),
#                   CallTaxaTableFromTSV('TS2010CLZM02.merged.tsv',level,"TS2010CLZM02"),
#                   CallTaxaTableFromTSV('TS2010CLZM03.merged.tsv',level,"TS2010CLZM03"),
#                   CallTaxaTableFromTSV('TS2010CLZM04.merged.tsv',level,"TS2010CLZM04"),
#                   CallTaxaTableFromTSV('TS2010CLZM05.merged.tsv',level,"TS2010CLZM05"),
#                   CallTaxaTableFromTSV('TS2010CLZM06.merged.tsv',level,"TS2010CLZM06"),
#                   CallTaxaTableFromTSV('TS2010CLZM07.merged.tsv',level,"TS2010CLZM07"),
#                   CallTaxaTableFromTSV('TS2010CLZM08.merged.tsv',level,"TS2010CLZM08"),
#                   CallTaxaTableFromTSV('TS2010CLZM09.merged.tsv',level,"TS2010CLZM09"),
#                   CallTaxaTableFromTSV('TS2010CLZM10.merged.tsv',level,"TS2010CLZM10"),
#                   CallTaxaTableFromTSV('TS2010CLZM11.merged.tsv',level,"TS2010CLZM11"),
#                   CallTaxaTableFromTSV('TS2010CLZM12.merged.tsv',level,"TS2010CLZM12"),
#                   CallTaxaTableFromTSV('TS2010CLZM13.merged.tsv',level,"TS2010CLZM13"),
                   CallTaxaTableFromTSV('TS2010CLZM14.merged.tsv',level,"TS2010CLZM14")))
#                   CallTaxaTableFromTSV('TS2010CLZM15.merged.tsv',level,"TS2010CLZM15"),
#                   CallTaxaTableFromTSV('TS2010CLZM16.merged.tsv',level,"TS2010CLZM16")))
  allsample[allsample$sample %in% c("TS2010CLZM01","TS2010CLZM02","TS2010CLZM03","TS2010CLZM14"),"group"]="bulk_SI_C"
  allsample[allsample$sample %in% c("TS2010CLZM04","TS2010CLZM05","TS2010CLZM06","TS2010CLZM15"),"group"]="bulk_LI_AB"
  allsample[allsample$sample %in% c("TS2010CLZM07","TS2010CLZM08","TS2010CLZM09","TS2010CLZM13"),"group"]="bulk_LI_C"
  allsample[allsample$sample %in% c("TS2010CLZM010","TS2010CLZM011","TS2010CLZM012","TS2010CLZM16"),"group"]="bulk_SI_AB"
  return(allsample)
}

test<-getBulkTaxaData("Family")
test<-group_by(test,group,taxa)
test<-dplyr::summarise(test,mean_count=mean(count))
bulkA<-dplyr::filter(test,group=="bulk_SI_C")[,2:3]
rownames(bulkA)<-bulkA$taxa

SMTA<-CreateSeuratFromSMT('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT01/ST2008PTPT01/outs',taxa.level = "family")
SMTA<-SMTA@assays$Phyloseq_level@counts %>% rowSums()
names(SMTA)<-names(SMTA) %>% gsub('-','_',.)
SMTA<-SMTA[SMTA>=10]
sharedtaxa<-intersect(bulkA$taxa,names(SMTA))

cortable<-bulkA[sharedtaxa,]
colnames(cortable)[2]<-'bulk'
cortable$SMT<-SMTA[sharedtaxa] %>% as.numeric()
cortable$SMT<-cortable$SMT %>% log()
cortable$bulk<-cortable$bulk %>% log()
microbe<-coplot("bulk","SMT",cortable)+ylab("SMT pseudo-bulk")+xlab("bulk RNA")+theme(text=element_text(size=20),plot.margin=margin(0.25,0.75,0.25,0.25,unit="cm"))

SMTA_host<-Load10X_Spatial('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT01/ST2008PTPT01/outs')
bulk_host<-read.csv('/data/LyuLin/Spatial/Workdirectory/BulkRNA_with_ST_expression-2020-11-23/mouse_bulk_rawcounts_matrix.csv',header=T,row.names=1)
colnames(bulk_host)<-c("SI_C_bulk1","SI_C_bulk2","SI_C_bulk3",
                        "LI_AB_bulk1","LI_AB_bulk2","LI_AB_bulk3",
                        "LI_C_bulk1","LI_C_bulk2","LI_C_bulk3",
                        "SI_AB_bulk1","SI_AB_bulk2","SI_AB_bulk3",
                        "LI_C_bulk_mix","SI_C_bulk_mix","LI_AB_bulk_mix","SI_AB_bulk_mix")
glens=getlength(rownames(bulk_host),"mm9","ensGene")
glens[is.na(glens)]<-glens[!is.na(glens)] %>% mean
bulk_host<-(bulk_host*1000/glens) %>% round()
bulk_host<-tibble::rownames_to_column(bulk_host,var = 'geneID')

readmapping<-read.delim('/data/LyuLin/Spatial/Workdirectory/BulkRNA_with_ST_expression-2020-11-23/gene_id_type_name_ensembl')
rownames(readmapping)<-readmapping$X
readmapping$X<-NULL
geneidlist<-bulk_host$geneID
genes<-readmapping[geneidlist,]
genes$gene_type<-NULL
genes<-tibble::rownames_to_column(genes,var = 'geneID')
colnames(genes)<-c("geneID","geneName")
bulk_host<-left_join(bulk_host,genes,by='geneID')
bulk_host<-bulk_host[,c("SI_C_bulk_mix","geneName")]
bulk_host<-group_by(bulk_host,geneName)
bulk_host<-dplyr::summarise(bulk_host,count=sum(SI_C_bulk_mix))
bulk_host<-bulk_host %>% as.data.frame()
rownames(bulk_host)<-bulk_host$geneName
sharedgenes<-intersect(bulk_host$geneName,rownames(SMTA_host))
cortable2<-data.frame(gene_name=sharedgenes,bulk=bulk_host[sharedgenes,"count"],SMT=rowSums(SMTA_host)[sharedgenes]%>%as.numeric())
cortable2<-dplyr::filter(cortable2,SMT>=3,bulk>=3)
cortable2$bulk<-log(cortable2$bulk)
cortable2$SMT<-log(cortable2$SMT)
host<-coplot("bulk","SMT",cortable2,alpha = 0.25)+ylab("SMT pseudo-bulk")+xlab("bulk RNA")+theme(text=element_text(size=20),plot.margin=margin(0.25,0.75,0.25,0.25,unit="cm"))
ggarrange(plotlist = list(microbe,host),nrow=2)
ggsave('/data/LyuLin/Spatial/Workdirectory/SMT_figures/Fig1C.pdf',width=4.45,height=8.83)
ggsave('/data/LyuLin/Spatial/Workdirectory/SMT_figures/FigS1.pdf')
