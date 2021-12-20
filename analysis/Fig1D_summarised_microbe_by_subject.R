source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')
source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')

level<-"species"
bulk_level<-"taxid1"

test<-CreateSeuratFromSMT(SMT.out.path="/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT01/ST2008PTPT01/outs",taxa.level = level)
test<-getPartition(test,partition.num=3,mode = "SMT",labels = c("up","middle","down"))
SI_upper<-data.frame(SI_upper=subset(test,SMTpartition=="up")@assays$Spatial@counts %>% rowSums())
SI_middle<-data.frame(SI_middle=subset(test,SMTpartition=="middle")@assays$Spatial@counts %>% rowSums())
SI_buttom<-data.frame(SI_buttom=subset(test,SMTpartition=="down")@assays$Spatial@counts %>% rowSums())

test2<-CreateSeuratFromSMT(SMT.out.path="/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT03/ST2008PTPT03/outs",taxa.level = level)
test2<-getPartition(test2,partition.num = 4,mode = "SMT",eps = 1.6,labels = c("debris","middle","down","up"))
LI_upper<-data.frame(LI_upper=subset(test2,SMTpartition=="up")@assays$Spatial@counts %>% rowSums())
LI_middle<-data.frame(LI_middle=subset(test2,SMTpartition=="middle")@assays$Spatial@counts %>% rowSums())
LI_buttom<-data.frame(LI_buttom=subset(test2,SMTpartition=="down")@assays$Spatial@counts %>% rowSums())

test3<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR1X/outs',taxa.level = level)
test3<-getPartition(test3,partition.num = 2,mode = "SMT",labels = c("tumor","normal"))
P53_normal1<-data.frame(P53_normal1=subset(test3,SMTpartition=="normal")@assays$Spatial@counts %>% rowSums())
P53_tumor1<-data.frame(P53_tumor1=subset(test3,SMTpartition=="tumor")@assays$Spatial@counts %>% rowSums())

test4<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR2X/outs',taxa.level = level)
test4<-getPartition(test4,partition.num = 2,mode = "SMT",labels = c("tumor","normal"))
P53_normal2<-data.frame(P53_normal2=subset(test4,SMTpartition=="normal")@assays$Spatial@counts %>% rowSums())
P53_tumor2<-data.frame(P53_tumor2=subset(test4,SMTpartition=="tumor")@assays$Spatial@counts %>% rowSums())

test5<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR3X/outs',taxa.level = level)
test5<-getPartition(test5,partition.num = 2,mode = "SMT",labels = c("tumor","normal"))
P52_normal1<-data.frame(P52_normal1=subset(test5,SMTpartition=="normal")@assays$Spatial@counts %>% rowSums())
P52_tumor1<-data.frame(P52_tumor1=subset(test5,SMTpartition=="tumor")@assays$Spatial@counts %>% rowSums())

test6<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR4X/outs',taxa.level = level)
test6<-getPartition(test6,partition.num = 2,mode = "SMT",labels = c("tumor","normal"))
P52_normal2<-data.frame(P52_normal2=subset(test6,SMTpartition=="normal")@assays$Spatial@counts %>% rowSums())
P52_tumor2<-data.frame(P52_tumor2=subset(test6,SMTpartition=="tumor")@assays$Spatial@counts %>% rowSums())

test7<-CallTaxaTableFromTSV('/data/LyuLin/Spatial/Workdirectory/BulkRNA-2020-11-12/blastraw/TS2010CLZM01.merged.tsv',taxa.level=bulk_level,scale.factor=0,sample.id="TS2010CLZM01")
test7[test7$taxa==paste0(level,"_"),"taxa"]<-paste0(level,"_Unknown")
rownames(test7)<-test7$taxa
test7<-test7[1]
rownames(test7)<-rownames(test7) %>% gsub("_","-",.)
colnames(test7)<-"SI_bulk1"

test8<-CallTaxaTableFromTSV('/data/LyuLin/Spatial/Workdirectory/BulkRNA-2020-11-12/blastraw/TS2010CLZM02.merged.tsv',taxa.level=bulk_level,scale.factor=0,sample.id="TS2010CLZM02")
test8[test8$taxa==paste0(level,"_"),"taxa"]<-paste0(level,"_Unknown")
rownames(test8)<-test8$taxa
test8<-test8[1]
rownames(test8)<-rownames(test8) %>% gsub("_","-",.)
colnames(test8)<-"SI_bulk2"

test9<-CallTaxaTableFromTSV('/data/LyuLin/Spatial/Workdirectory/BulkRNA-2020-11-12/blastraw/TS2010CLZM03.merged.tsv',taxa.level=bulk_level,scale.factor=0,sample.id="TS2010CLZM03")
test9[test9$taxa==paste0(level,"_"),"taxa"]<-paste0(level,"_Unknown")
rownames(test9)<-test9$taxa
test9<-test9[1]
rownames(test9)<-rownames(test9) %>% gsub("_","-",.)
colnames(test9)<-"SI_bulk3"

test10<-CallTaxaTableFromTSV('/data/LyuLin/Spatial/Workdirectory/BulkRNA-2020-11-12/blastraw/TS2010CLZM07.merged.tsv',taxa.level=bulk_level,scale.factor=0,sample.id="TS2010CLZM07")
test10[test10$taxa==paste0(level,"_"),"taxa"]<-paste0(level,"_Unknown")
rownames(test10)<-test10$taxa
test10<-test10[1]
rownames(test10)<-rownames(test10) %>% gsub("_","-",.)
colnames(test10)<-"LI_bulk1"

test11<-CallTaxaTableFromTSV('/data/LyuLin/Spatial/Workdirectory/BulkRNA-2020-11-12/blastraw/TS2010CLZM08.merged.tsv',taxa.level=bulk_level,scale.factor=0,sample.id="TS2010CLZM08")
test11[test11$taxa==paste0(level,"_"),"taxa"]<-paste0(level,"_Unknown")
rownames(test11)<-test11$taxa
test11<-test11[1]
rownames(test11)<-rownames(test11) %>% gsub("_","-",.)
colnames(test11)<-"LI_bulk2"

test12<-CallTaxaTableFromTSV('/data/LyuLin/Spatial/Workdirectory/BulkRNA-2020-11-12/blastraw/TS2010CLZM09.merged.tsv',taxa.level=bulk_level,scale.factor=0,sample.id="TS2010CLZM09")
test12[test12$taxa==paste0(level,"_"),"taxa"]<-paste0(level,"_Unknown")
rownames(test12)<-test12$taxa
test12<-test12[1]
rownames(test12)<-rownames(test12) %>% gsub("_","-",.)
colnames(test12)<-"LI_bulk3"

df<-Reduce(mergeByRownames,list(SI_upper,SI_middle,SI_buttom,LI_upper,LI_middle,LI_buttom,P53_normal1,P53_normal2,P53_tumor1,P53_tumor2,P52_normal1,P52_normal2,P52_tumor1,P52_tumor2,test7,test8,test9,test10,test11,test12))
df[is.na(df)]<-0
rownames(df)<-rownames(df) %>% gsub("^.*-","",.) %>% gsub("^.*_","",.)
taxadf<-read.delim('./taxa.level.matrix.tsv',row.names=1,header = F)
taxadf<-rownames_to_column(taxadf,var="taxid")
colnames(taxadf)[2:9]<-c("superkingdom","kingdom","phylum","class","order","family","genus","species")
df<-rownames_to_column(df,var="taxid")
df<-left_join(df,taxadf,by="taxid")
saveRDS(df,'/data/LyuLin/Spatial/Workdirectory/SMT_figures/Fig1D_taxid_data_for_RFeng.rds')
saveRDS(df,'/data/LyuLin/Spatial/Workdirectory/SMT_figures/Fig1D_order_data_for_RFeng.rds')

df<-readRDS('/data/LyuLin/Spatial/Workdirectory/SMT_figures/Fig1D_data_for_RFeng.rds')
df<-TransformDataframe(df,0.05)
df$features<-df$features %>% gsub("^.*-","",.) %>% gsub("^.*_","",.)
ggplot(df)+geom_bar(aes(x=spot,y=values,fill=features),stat="identity",position="fill")+scale_fill_d3("category20")+theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1))
