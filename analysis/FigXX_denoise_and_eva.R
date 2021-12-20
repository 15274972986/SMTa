source('/data/LyuLin/Scripts/spatial_scripts/Spatial_core_functions.R')
source('/data/LyuLin/Scripts/spatial_scripts/SMT_core_functions.R')
setwd('/data/LyuLin/Spatial/Workdirectory/SMT_figures/')

mmSI<-CreateSeuratFromSMT('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT01/ST2008PTPT01/outs',taxa.level="species")
mmLI<-CreateSeuratFromSMT('/data/LyuLin/Spatial/Workdirectory/2020-8-24-Spatial_data_and_run/2020-8-24-Spatial_ST2008PTPT03/ST2008PTPT03/outs',taxa.level="species")

mmSISpotSig0<-CalSpotSignificance(mmSI,n.permutation=1000)
mmSISpotSig<-mmSISpotSig0 %>% dplyr::filter(.,p.value<0.01)
mmSIUMISig0<-CalUMISignificance(mmSI,n.permutation=1000)
mmSIUMISig<-mmSIUMISig0 %>% dplyr::filter(.,p.value<0.01)
valid.taxa.mmSI<-unique(mmSISpotSig$taxa,mmSIUMISig$taxa)
valid.taxa.mmSI.raw<-FetchData(mmSI,vars=valid.taxa.mmSI)
valid.taxa.mmSI.spots<-rownames(valid.taxa.mmSI.raw[rowSums(valid.taxa.mmSI.raw)>=1,])
invalid.taxa.mmSI<-mmSISpotSig0$taxa[!(mmSISpotSig0$taxa %in% valid.taxa.mmSI)]
invalid.taxa.mmSI.raw<-FetchData(mmSI,vars=invalid.taxa.mmSI)
invalid.taxa.mmSI.spots<-rownames(invalid.taxa.mmSI.raw[rowSums(invalid.taxa.mmSI.raw)>=1,])
sharedspots<-intersect(valid.taxa.mmSI.spots,invalid.taxa.mmSI.spots)
valid.taxa.mmSI.spots<-setdiff(valid.taxa.mmSI.spots,sharedspots)
invalid.taxa.mmSI.spots<-setdiff(invalid.taxa.mmSI.spots,sharedspots)
mmSI<-getPeriTissueLayers(mmSI,Cells(mmSI))
valid.taxa.mmSI.meta<-mmSI@meta.data[valid.taxa.mmSI.spots,]
valid.taxa.mmSI.meta$valididty<-"valid taxa"
valid.taxa.mmSI.meta<-rownames_to_column(valid.taxa.mmSI.meta,"barcode")
invalid.taxa.mmSI.meta<-mmSI@meta.data[invalid.taxa.mmSI.spots,]
invalid.taxa.mmSI.meta$valididty<-"invalid taxa"
invalid.taxa.mmSI.meta<-rownames_to_column(invalid.taxa.mmSI.meta,"barcode")
base.mmSI.meta<-mmSI@meta.data
base.mmSI.meta$valididty<-"base"
base.mmSI.meta<-rownames_to_column(base.mmSI.meta,"barcode")
data<-base::Reduce(rbind,list(valid.taxa.mmSI.meta,invalid.taxa.mmSI.meta,base.mmSI.meta))


ggplot(data)+geom_boxplot(aes(x=valididty,y=dist.from.tissue),width=0.6,outlier.alpha=0)+
  geom_signif(aes(x=valididty,y=dist.from.tissue),step_increase = 0.1,map_signif_level = F,
              comparisons=list(c("valid taxa","invalid taxa"),
                               c("valid taxa","base"),
                               c("base","invalid taxa")),test ="wilcox.test")+
  ylab("distance to tissue")+
  theme_cowplot(font_size = 15)

SMTPlotMod(mmSI,"species-Rhodotorula-glutinis",alpha = c(1,1),colors=c("orange","red"),title="Rhodotorula glutinis")
SMTPlotMod(mmSI,"species-Diversispora-versiformis",alpha = c(1,1),colors=c("orange","red"),title="Diversispora versiformis")
SMTPlotMod(mmSI,"species-Inocybe-dulcamara",alpha = c(1,1),colors=c("orange","red"),title="Inocybe dulcamara")

data.frame(species=)

#method compare
venn.data<-list("spot_based"=mmSISpotSig$taxa,"UMI_based"=mmSIUMISig$taxa)
plot_grid(venn.diagram(venn.data,fill=pal_d3("category10")(2),filename = NULL))+coord_fixed()

#mmLI denoise
mmLISpotSig0<-CalSpotSignificance(mmLI,n.permutation=1000)
mmLISpotSig<-mmLISpotSig0 %>% dplyr::filter(.,p.value<0.01)
mmLIUMISig0<-CalUMISignificance(mmLI,n.permutation=1000)
mmLIUMISig<-mmLIUMISig0 %>% dplyr::filter(.,p.value<0.01)
valid.taxa.mmLI<-unique(mmLISpotSig$taxa,mmLIUMISig$taxa)
valid.taxa.mmLI.raw<-FetchData(mmLI,vars=valid.taxa.mmLI)
valid.taxa.mmLI.spots<-rownames(valid.taxa.mmLI.raw[rowSums(valid.taxa.mmLI.raw)>=1,])
invalid.taxa.mmLI<-mmLISpotSig0$taxa[!(mmLISpotSig0$taxa %in% valid.taxa.mmLI)]
invalid.taxa.mmLI.raw<-FetchData(mmLI,vars=invalid.taxa.mmLI)
invalid.taxa.mmLI.spots<-rownames(invalid.taxa.mmLI.raw[rowSums(invalid.taxa.mmLI.raw)>=1,])
sharedspots2<-intersect(valid.taxa.mmLI.spots,invalid.taxa.mmLI.spots)
valid.taxa.mmLI.spots<-setdiff(valid.taxa.mmLI.spots,sharedspots2)
invalid.taxa.mmLI.spots<-setdiff(invalid.taxa.mmLI.spots,sharedspots2)
mmLI<-getPeriTissueLayers(mmLI,Cells(mmLI))
valid.taxa.mmLI.meta<-mmLI@meta.data[valid.taxa.mmLI.spots,]
valid.taxa.mmLI.meta$valididty<-"valid taxa"
valid.taxa.mmLI.meta<-rownames_to_column(valid.taxa.mmLI.meta,"barcode")
invalid.taxa.mmLI.meta<-mmLI@meta.data[invalid.taxa.mmLI.spots,]
invalid.taxa.mmLI.meta$valididty<-"invalid taxa"
invalid.taxa.mmLI.meta<-rownames_to_column(invalid.taxa.mmLI.meta,"barcode")
base.mmLI.meta<-mmLI@meta.data
base.mmLI.meta$valididty<-"base"
base.mmLI.meta<-rownames_to_column(base.mmLI.meta,"barcode")
data2<-base::Reduce(rbind,list(valid.taxa.mmLI.meta,invalid.taxa.mmLI.meta,base.mmLI.meta))
ggplot(data2)+geom_boxplot(aes(x=valididty,y=dist.from.tissue),width=0.6,outlier.alpha=0)+
  geom_signif(aes(x=valididty,y=dist.from.tissue),step_increase = 0.1,map_signif_level = F,
              comparisons=list(c("valid taxa","invalid taxa"),
                               c("valid taxa","base"),
                               c("base","invalid taxa")),test ="wilcox.test")+
  ylab("distance to tissue")+
  theme_cowplot(font_size = 15)

SMTPlotMod(mmLI,"species-Faecalibaculum-rodentium",alpha = c(1,1),colors=c("orange","red"),title=NULL)

hsP53A<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR1X/outs',taxa.level="species")
hsP53B<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR2X/outs',taxa.level="species")
hsP52C<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR3X/outs',taxa.level="species")
hsP52D<-CreateSeuratFromSMT('/data/LyuLin/Spatial/fastqs/Spatial_2nd_trial/ST2103ZJCR4X/outs',taxa.level="species")
hsP53ASpotSig0<-CalSpotSignificance(hsP53A,n.permutation=1000)
hsP53ASpotSig<-hsP53ASpotSig0 %>% dplyr::filter(.,p.value<0.01)
hsP53AUMISig0<-CalUMISignificance(hsP53A,n.permutation=1000)
hsP53AUMISig<-hsP53AUMISig0 %>% dplyr::filter(.,p.value<0.01)
valid.taxa.hsP53A<-unique(hsP53ASpotSig$taxa,hsP53AUMISig$taxa)

hsP53BSpotSig0<-CalSpotSignificance(hsP53B,n.permutation=1000)
hsP53BSpotSig<-hsP53BSpotSig0 %>% dplyr::filter(.,p.value<0.01)
hsP53BUMISig0<-CalUMISignificance(hsP53B,n.permutation=1000)
hsP53BUMISig<-hsP53BUMISig0 %>% dplyr::filter(.,p.value<0.01)
valid.taxa.hsP53B<-unique(hsP53BSpotSig$taxa,hsP53BUMISig$taxa)

hsP52CSpotSig0<-CalSpotSignificance(hsP52C,n.permutation=1000)
hsP52CSpotSig<-hsP52CSpotSig0 %>% dplyr::filter(.,p.value<0.01)
hsP52CUMISig0<-CalUMISignificance(hsP52C,n.permutation=1000)
hsP52CUMISig<-hsP52CUMISig0 %>% dplyr::filter(.,p.value<0.01)
valid.taxa.hsP52C<-unique(hsP52CSpotSig$taxa,hsP52CUMISig$taxa)

hsP52DSpotSig0<-CalSpotSignificance(hsP52D,n.permutation=1000)
hsP52DSpotSig<-hsP52DSpotSig0 %>% dplyr::filter(.,p.value<0.01)
hsP52DUMISig0<-CalUMISignificance(hsP52D,n.permutation=1000)
hsP52DUMISig<-hsP52DUMISig0 %>% dplyr::filter(.,p.value<0.01)
valid.taxa.hsP52D<-unique(hsP52DSpotSig$taxa,hsP52DUMISig$taxa)

hsP70B<-CreateSeuratFromSMT('/data/LyuLin/Scripts/spacewrapper3.0/output/ST2110CLLX2X/ST2110CLLX2X/outs/',taxa.level="species")
hsP70BSpotSig0<-CalSpotSignificance(hsP70B,n.permutation=1000)
hsP70BSpotSig<-hsP70BSpotSig0 %>% dplyr::filter(.,p.value<0.01)
hsP70BUMISig0<-CalUMISignificance(hsP70B,n.permutation=1000)
hsP70BUMISig<-hsP70BUMISig0 %>% dplyr::filter(.,p.value<0.01)
valid.taxa.hsP70B<-unique(hsP70BSpotSig$taxa,hsP70BUMISig$taxa)

hsP70C<-CreateSeuratFromSMT('/data/LyuLin/Scripts/spacewrapper3.0/output/ST2110CLLX3X/ST2110CLLX3X/outs/',taxa.level="species")
hsP70CSpotSig0<-CalSpotSignificance(hsP70C,n.permutation=1000)
hsP70CSpotSig<-hsP70CSpotSig0 %>% dplyr::filter(.,p.value<0.01)
hsP70CUMISig0<-CalUMISignificance(hsP70C,n.permutation=1000)
hsP70CUMISig<-hsP70CUMISig0 %>% dplyr::filter(.,p.value<0.01)
valid.taxa.hsP70C<-unique(hsP70CSpotSig$taxa,hsP70CUMISig$taxa)

hsP70D<-CreateSeuratFromSMT('/data/LyuLin/Scripts/spacewrapper3.0/output/ST2110CLLX4X/ST2110CLLX4X/outs/',taxa.level="species")
hsP70DSpotSig0<-CalSpotSignificance(hsP70D,n.permutation=1000)
hsP70DSpotSig<-hsP70DSpotSig0 %>% dplyr::filter(.,p.value<0.01)
hsP70DUMISig0<-CalUMISignificance(hsP70D,n.permutation=1000)
hsP70DUMISig<-hsP70DUMISig0 %>% dplyr::filter(.,p.value<0.01)
valid.taxa.hsP70D<-unique(hsP70DSpotSig$taxa,hsP70DUMISig$taxa)

hsP70DSpotSig0 %>% dplyr::filter(.,taxa=="species-Schizosaccharomyces-pombe")
SMTPlot(hsP52C,"nCount_Spatial",alpha=0)+NoLegend()+SMTPlotMod(hsP52C,"species-Schizosaccharomyces-pombe",alpha = c(1,1),title = NULL)
