CellDimPlot(
  srt = bladder_merge_afterdouble_annotation, group.by = "celltype",
  reduction = "UMAP", theme_use = "theme_blank",palette = "Paired"
)
CellDimPlot(
  srt = bladder_merge_afterdouble_annotation, group.by = "celltype", stat.by = "Phase",
  reduction = "UMAP", theme_use = "theme_blank"
)
FeatureDimPlot(
  srt = bladder_merge_afterdouble_annotation, features = c("Epcam","Krt19","Col1a2", "Ptprc"),
  reduction = "UMAP", theme_use = "theme_blank"
)
ht <- GroupHeatmap(
  srt = bladder_merge_afterdouble_annotation,
  features = top5$gene,
  group.by = c("celltype"),
  heatmap_palette = "YlOrRd",
  #cell_annotation = c("Phase", "G2M_score", "Cdh2"),
  cell_annotation_palette = "Set2",
  show_row_names = TRUE, row_names_side = "left",
  add_dot = TRUE, add_reticle = TRUE
)
print(ht$plot)



bladder_merge_afterdouble_annotation <- RunDEtest(srt = bladder_merge_afterdouble_annotation, group_by = "celltype", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = bladder_merge_afterdouble_annotation, group_by = "orig.ident")
DEGs <- bladder_merge_afterdouble_annotation@tools$DEtest_celltype$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 0.25 & p_val_adj < 0.05), ]
# Annotate features with transcription factors and surface proteins
bladder_merge_afterdouble_annotation <- AnnotateFeatures(bladder_merge_afterdouble_annotation, species = "Mus_musculus", db = c("TF", "CSPA"))
ht <- FeatureHeatmap(
  srt = bladder_merge_afterdouble_annotation, group.by = "celltype", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Mus_musculus", db = c("GO_BP"), anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4,topTerm = 3,heatmap_palette = "RdBu",cell_annotation_palette =  brewer.pal(6,"Paired"), anno_features = F
)
print(ht$plot)
bladder_merge_afterdouble_annotation <- RunEnrichment(
  srt = bladder_merge_afterdouble_annotation, group_by = "orig.ident", db=c("GO_BP","GO_CC","GO_MF"), species = "Mus_musculus",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
)
EnrichmentPlot(
  srt = bladder_merge_afterdouble_annotation, group_by = "orig.ident",
  plot_type = "lollipop",topTerm = 4,db=c("GO_BP","GO_CC","GO_MF")
)
bladder_merge_afterdouble_annotation <- RunSlingshot(srt = bladder_merge_afterdouble_annotation, group.by = "celltype_sub", reduction = "UMAP",start = "basal_1")
FeatureDimPlot(bladder_merge_afterdouble_annotation, features = paste0("Lineage", 1:2), reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(bladder_merge_afterdouble_annotation, group.by = "seurat_clusters", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1,split.by = "orig.ident")

FeatureDimPlot(
  srt = bladder_merge_afterdouble_annotation, features = c( "Malat1","Nkx3-1"),
  compare_features = T, label = T, label_insitu = TRUE,
  reduction = "umap", theme_use = "theme_blank"
)
bladder_merge_afterdouble_annotation <- RunDynamicFeatures(srt = bladder_merge_afterdouble_annotation ,lineages = c("Lineage1", "Lineage2"), n_candidates = 200)
Dynamic <- DynamicHeatmap(
  srt = bladder_merge_afterdouble_annotation, lineages = c("Lineage1","Lineage2"),
  use_fitted = TRUE, n_split = 2, reverse_ht = "Lineage1",
  species = "Mus_musculus", db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  heatmap_palette = "viridis", cell_annotation = "orig.ident",
  separate_annotation = list("celltype_sub",c("Egf","Fos")), separate_annotation_palette = c("Paired", "Set1"),
  feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
)
print(Dynamic$plot)

bladder_merge_afterdouble_annotation <- Standard_SCP(bladder_merge_afterdouble_annotation)
CellDimPlot(bladder_merge_afterdouble_annotation, group.by = c("Standardclusters", "celltype_sub"), label = TRUE, theme_use = "theme_blank")
bladder_merge_afterdouble_annotation <- RunMonocle3(srt = bladder_merge_afterdouble_annotation, annotation = "celltype_sub", clusters = "Standardclusters")
trajectory <- bladder_merge_afterdouble_annotation@tools$Monocle3$trajectory
p1 <- CellDimPlot(bladder_merge_afterdouble_annotation, group.by = "Monocle3_partitions", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
p2 <- CellDimPlot(bladder_merge_afterdouble_annotation, group.by = "Monocle3_clusters", reduction = "StandardUMAP2D", label = TRUE, theme_use = "theme_blank") + trajectory
p3 <- FeatureDimPlot(bladder_merge_afterdouble_annotation, features = "Monocle3_Pseudotime", reduction = "StandardUMAP2D", theme_use = "theme_blank") + trajectory
print(p1 + p2 + p3)

bladder_merge_afterdouble_annotation <- RunMonocle2(srt = bladder_merge_afterdouble_annotation, annotation = "celltype_sub")
names(bladder_merge_afterdouble_annotation@tools$Monocle2)
trajectory <- bladder_merge_afterdouble_annotation@tools$Monocle2$trajectory

p1 <- CellDimPlot(bladder_merge_afterdouble_annotation, group.by = "Monocle2_State", reduction = "DDRTree", label = TRUE, theme_use = "theme_blank") + trajectory
CellDimPlot(bladder_merge_afterdouble_annotation, group.by = "Monocle2_State", reduction = "UMAP", label = TRUE, theme_use = "theme_blank",palette = "Spectral")
p3 <- FeatureDimPlot(bladder_merge_afterdouble_annotation, features = "Monocle2_Pseudotime", reduction = "UMAP", theme_use = "theme_blank")
print(p1 + p2 + p3)

CellDimPlot3D(srt = prostate_merge_abtainHb_afterdoulet, group.by = "orig.ident")
EnrichmentPlot(
  srt = bladder_merge_afterdouble_annotation, group_by = "orig.ident",
  plot_type = "bar"
)
EnrichmentPlot(
  srt = bladder_merge_afterdouble_annotation, group_by = "orig.ident",
  plot_type = "wordcloud"
)
EnrichmentPlot(
  srt = bladder_merge_afterdouble_annotation, group_by = "orig.ident",
  plot_type = "wordcloud", word_type = "feature"
)
EnrichmentPlot(
  srt = bladder_merge_afterdouble_annotation, group_by = "orig.ident",
  plot_type = "network"
)
EnrichmentPlot(
  srt = bladder_merge_afterdouble_annotation, group_by = "orig.ident",
  plot_type = "enrichmap"
)
