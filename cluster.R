prostate_merge_afterdoublefinder <- NormalizeData(prostate_merge_afterdoublefinder)
prostate_merge_afterdoublefinder <- FindVariableFeatures(prostate_merge_afterdoublefinder, selection.method = "vst", nfeatures = 2000)
#plotgene1 <- VariableFeaturePlot(prostate_merge_afterdoublefinder)
#plotgene2 <- LabelPoints(plot = plotgene1, points = top10, repel = TRUE)
#prostate_merge_afterdoublefinder <- SCTransform(prostate_merge_afterdoublefinder, verbose = T, vars.to.regress = c("nCount_RNA", "percent.mt"), conserve.memory = T)


prostate_merge_afterdoublefinder <- ScaleData(prostate_merge_afterdoublefinder)
prostate_merge_afterdoublefinder <- RunPCA(prostate_merge_afterdoublefinder, features = VariableFeatures(prostate_merge_afterdoublefinder),npcs = 100)

library(harmony)
prostate_merge_afterdoublefinder <- RunHarmony(prostate_merge_afterdoublefinder, "batch")

ElbowPlot(prostate_merge_afterdoublefinder,ndims = 100)
prostate_merge_afterdoublefinder <- JackStraw(prostate_merge_afterdoublefinder, num.replicate = 100)
prostate_merge_afterdoublefinder <- ScoreJackStraw(prostate_merge_afterdoublefinder, dims = 1:100)
JackStrawPlot(prostate_merge_afterdoublefinder, dims = 1:100)
prostate_merge_afterdoublefinder <- FindNeighbors(prostate_merge_afterdoublefinder, reduction="pca",dims = 1:30)
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
prostate_merge_afterdoublefinder <- RunUMAP(prostate_merge_afterdoublefinder,reduction = "pca", dims = 1:30)
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






