prostate_merge_afterdoublefinder_afterdoublefinder <-  merge(Aged_AP_1_afterdoublefinder, 
                                           y = c(Aged_AP_2_afterdoublefinder,Aged_VDLP_1_afterdoublefinder,Aged_VDLP_2_afterdoublefinder,Young_AP_1_afterdoublefinder,Young_AP_2_afterdoublefinder,Young_VDLP_1_afterdoublefinder,Young_VDLP_2_afterdoublefinder), add.cell.ids = c("1", "2","3","4","5","6","7","8"), project = "big",merge.data=TRUE)
sc_data <- SweepFeatureData(object = prostate_merge_afterdoublefinder, features = c("Hbb-bs","Hba-a2","Hba-a1","Hbb-bt"))
which(rownames(counts) == "Hbb-bs") | which(rownames(counts) == "Hba-a2")| which(rownames(counts) == "Hba-a1") |which(rownames(counts) == "Hbb-bt")
prostate_merge_afterdoublefinder
prostate_merge_afterdoublefinder <- NormalizeData(prostate_merge_afterdoublefinder)
prostate_merge_afterdoublefinder <- FindVariableFeatures(prostate_merge_afterdoublefinder, selection.method = "vst", nfeatures = 2000)
#plotgene1 <- VariableFeaturePlot(prostate_merge_afterdoublefinder)
#plotgene2 <- LabelPoints(plot = plotgene1, points = top10, repel = TRUE)
#prostate_merge_afterdoublefinder <- SCTransform(prostate_merge_afterdoublefinder, verbose = T, vars.to.regress = c("nCount_RNA", "percent.mt"), conserve.memory = T)


prostate_merge_afterdoublefinder <- ScaleData(prostate_merge_afterdoublefinder)
prostate_merge_afterdoublefinder <- RunPCA(prostate_merge_afterdoublefinder, features = VariableFeatures(prostate_merge_afterdoublefinder),npcs = 100)

library(harmony)
prostate_merge_afterdoublefinder <- RunHarmony(prostate_merge_afterdoublefinder, "orig.ident")

ElbowPlot(prostate_merge_afterdoublefinder,ndims = 100)
prostate_merge_afterdoublefinder <- JackStraw(prostate_merge_afterdoublefinder, num.replicate = 100)
prostate_merge_afterdoublefinder <- ScoreJackStraw(prostate_merge_afterdoublefinder, dims = 1:100)
JackStrawPlot(prostate_merge_afterdoublefinder, dims = 1:100)
prostate_merge_afterdoublefinder <- FindNeighbors(prostate_merge_afterdoublefinder, reduction="pca",dims = 1:50)
prostate_merge_afterdoublefinder <- FindNeighbors(prostate_merge_afterdoublefinder, reduction="harmony",dims = 1:30)
res.used <- seq(0.1,1,by=0.1)
res.used
for(i in res.used){
  prostate_merge_afterdoublefinder <- FindClusters(prostate_merge_afterdoublefinder,  resolution = res.used)
}



library(clustree)
clus.tree.out <- clustree(prostate_merge_afterdoublefinder) +
  theme(legend.position = "bottom") + 
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(low = "grey80", high = "red")
clus.tree.out

prostate_merge_afterdoublefinder <- FindClusters(prostate_merge_afterdoublefinder, resolution = 0.1)
prostate_merge_afterdoublefinder <- RunUMAP(prostate_merge_afterdoublefinder,reduction = "pca", dims = 1:50)
prostate_merge_afterdoublefinder <- RunUMAP(prostate_merge_afterdoublefinder,reduction = "harmony", dims = 1:30)
DimPlot(prostate_merge_afterdoublefinder, reduction = "umap",label = T)

pdf("afterharmony.pdf")
DimPlot(prostate_merge_afterdoublefinder, reduction = "umap",label = F,group.by = "orig.ident")
dev.off()

prostate_merge_afterdoublefinder <- RunTSNE(prostate_merge_afterdoublefinder,reduction = "harmony",dims = 1:30)
DimPlot(prostate_merge_afterdoublefinder, reduction = "tsne",label = T,pt.size = 1)
library(RColorBrewer)
colorcount <- length(unique(prostate_merge_afterdoublefinder@meta.data$seurat_clusters))
getPalette <- colorRampPalette(brewer.pal(12,"Paired"))
celltype_colors <- getPalette(colorcount)
DimPlot(prostate_merge_afterdoublefinder, reduction = "tsne",label = T,cols = celltype_colors)


cell_markers<- FindAllMarkers(prostate_merge_afterdoublefinder_3, 
                                              only.pos = TRUE,  # 只返回positive基因
                                              min.pct = 0.25,logfc.threshold = 0.25
) #只计算至少在(两簇细胞总数的)25%的细胞中有表达的基因
saveRDS(prostate_merge_afterdoublefinder, file = "prostate_merge_afterdoublefinder.RDS")
save(cell_markers,file = "cellmarkers.Rdata")


cell_markers_3 %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top30 <- cell_markers_3 %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

genes_of_interest <- top30$gene


avg_expr <- AverageExpression(object = prostate_merge_afterdoublefinder_3, genes.use = genes_of_interest)


DoHeatmap(object =prostate_merge_afterdoublefinder_3, features = avg_expr$gene, group.by = "celltype",label = F) + 
  scale_fill_viridis()

DoHeatmap(prostate_merge_afterdoublefinder_3, features = top10$gene) + NoLegend()

fibr <- c("S100A4","COL3A1","COL1A1","COL1A2","ACTA2","VIM",
          "DCN","FBLN1","LUM","PI16","NBL1","MFAP4")
epi <- c("KRT18","KRT19","KRT17","UPK1A","UPK1B","UPK3B",
         "UPK3A","Krt7","Krt15")
SMC <- c("TAGLN","DES","CNN1","ACTG2","TPM2","MYL9","MYH11","MYLK")
endo <- c("CALD1","RRAD","MUSTN1","SELE","PECAM1","VCAM1"  
          ,"CDH5","PLVAP","APQ1")
infla <- c("CD3E","CD79A")
NC <- c("GPM6A")
szcell <- c("KRT20","UPK1A","UPK1B","UPK3B","UPK3A","UPK2")
zjcell <- c("KRT13","UPK1A","UPK1B","UPK3B","UPK3A","TP63","SHH")
basal <- c("KRT5","KRT17","TP63","SHH")
"KRT1","KRT10","DSG1","CDH1"
"DSC1","KRT2","IVL","TGM3"
##0,1,5,11,14 epi
##2,10,12  fibr smc endo
##3  fibr endo
##4,8,9,13  fibr endo
##6,16  infla
##7,15,18,19 2,10,3,4,8,13,17 fibr
##17  fibr epi
##12 smc
##9 endo
clus <- 0:3
for(i in 1:length(clus)){print(table(merge.markers$gene[merge.markers$cluster==clus]%in% fibr))}
table(metrge.markers$gene[merge.markers$cluster==clus]%in% fibr)

FeaturePlot(blad_sub, features = "CA9",reduction = "tsne")
VlnPlot(prostate_merge_afterdoublefinder,features = "CA9",pt.size=0)
VlnPlot(prostate_merge_afterdoublefinder,features = epi,pt.size=0)
VlnPlot(prostate_merge_afterdoublefinder,features = SMC,pt.size=0)
DotPlot()

BiocManager::install("SingleR")
BiocManager::install("celldex")
library(SingleR)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
meta2=prostate_merge_afterdoublefinder@meta.data #prostate_merge_afterdoublefinder的meta文件，包含了seurat的聚类结果
head(meta2)
merge_for_SingleR <- GetAssayData(prostate_merge_afterdoublefinder, slot="data") ##获取标准化矩阵
merge.hesc <- SingleR(test = merge_for_SingleR, ref = hpca.se, labels = hpca.se$label.main) #
merge.hesc
#seurat 和 singleR的table表
table(merge.hesc$labels,meta2$seurat_clusters)
prostate_merge_afterdoublefinder@meta.data$labels <-merge.hesc$labels
print(DimPlot(prostate_merge_afterdoublefinder, group.by = c("seurat_clusters", "labels"),reduction = "umap",label = T))
top10 <- merge.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(prostate_merge_afterdoublefinder, features = top10$gene) + NoLegend()
FeaturePlot(prostate_merge_afterdoublefinder, features = c("S100A2","FABP4","KRT17" ))






