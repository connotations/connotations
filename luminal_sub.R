matrix <- subset(prostate_merge_exceptseminal, celltype=="luminal"|celltype=="luminal progenitor")
matrix <- NormalizeData(matrix)
matrix <- FindVariableFeatures(matrix, selection.method = "vst", nfeatures = 2000)
#scale.genes <-  rownames(matrix)
matrix <- ScaleData(matrix)
#matrix <- SCTransform(matrix, verbose = T, vars.to.regress = c("nCount_RNA", "percent.mt"), conserve.memory = T)
matrix <- RunPCA(matrix, features = VariableFeatures(matrix),npcs = 50)
matrix <- RunHarmony(matrix, "orig.ident")
ElbowPlot(matrix,ndims = 50)
matrix <- FindNeighbors(matrix,dims = 1:30)
matrix <- FindClusters(matrix,resolution = 0.8)
matrix <- RunUMAP(matrix,reduction = "pca",dims = 1:30)
matrix <- RunTSNE(matrix,reduction = "harmony",dims = 1:15)
DimPlot(matrix,reduction = "tsne",label = T,group.by = "celltype")
DimPlot(matrix,reduction = "umap",label = T,group.by = "celltype")
DimPlot(matrix,reduction = "umap",label = F,group.by = "orig.ident")
DimPlot(matrix,reduction = "umap",label = F,split.by = "group")

luminal_age_markers_2 <- FindMarkers(object = matrix, group.by = "orig.ident",ident.1="Old",ident.2="Young")
av <- AverageExpression(matrix,group.by="seurat_clusters",assay="RNA")
av=av[[1]]
cg <- names(tail(sort(apply(av,1,sd)),1000))
pheatmap::pheatmap(cor(av[cg,],method="spearman"))
matrix_markers <- FindAllMarkers(matrix, 
                                 only.pos = TRUE,  # 只返回positive基因
                                 min.pct = 0.25,logfc.threshold = 0.25)
top20 <- matrix_markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC)
FeaturePlot(matrix,c("Pate4","Pde4d","Ren1","Ly6a","Oasl2","Ly6e","Malat1","Alas2","Wfdc10","Oit1","Nrg1","Fxyd4","Msmb"))
DotPlot(matrix,features = c("PDCD1","LAG3","CTLA4","TIGIT"))
DotPlot(matrix,features = c("CCR7","SELL"))
DotPlot(matrix,features = c("MKI67","TOP2A"))
DotPlot(matrix,features = c("GZMK","CXCR4","CXCR3","CD44"))
DotPlot(matrix,features = c("CD6","XCL1","XCL2"))
DotPlot(matrix,features = c("KLRG1","CX3CR1","FCGR3A"))
DotPlot(matrix,features = c("IGKC"))
FeaturePlot(matrix,features = c("Mki67"))
geneset <- c("CD3D","CD4","CD8A",
             "HAVCR2","LAG3","PDCD1","CTLA4","TIGIT","BTLA","KLRC1",
             "CCR7","LEF1","SELL","TCF7",
             "GNLY","IFNG","NKG7","PRF1","GZMA","GZMB","GZMH","GZMK",
             
             "ANXA1","ANKRD28","IL7R","CD69","CD40LG",
             "FOXP3","IL2RA","IKZF2",
             "NCR1","NCAM1","TYROBP","FGFBP2","KLRD1","KLRF1","KLRB1","CX3CR1","FCGR3A")
#计算平均表达量
gene_cell_exp <- AverageExpression(matrix,
                                   features = geneset,
                                   group.by = 'seurat_clusters',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)

#complexheatmap作图
library(ComplexHeatmap)
#顶部细胞类型注释
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c('0'="#9ECABE",
                                                  '1'="#F6F5B4",
                                                  '2'="#2F528F",
                                                  "3"="#E3AD68",
                                                  "4"="#ACD45E",
                                                  "5"="",
                                                  "6"="",
                                                  "7"="",
                                                  "8"="",
                                                  "9"="",
                                                  "10"="",
                                                  "11"="",
                                                  "12"="",
                                                  "13"="",
                                                  "14"="",
                                                  "15"="",
                                                  "16"="",
                                                  "17"="",
                                                  "18"="")))#颜色设置
#数据标准化缩放一下
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
marker_exp <- marker_exp[,c(1,8,7,17,12,3,4,5,6,16,15,19,10,13,2,9,11,14,18)]
for (i in 1:39) {for(j in 1:19){if (marker_exp[i,j]>2){marker_exp[i,j]=2} }}
test <- marker_exp[,c(1,8,3,4,5,6,16,10,13,2,9,11,14,18)]
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = colorRampPalette(c("#9999FF","white","#AA0000"))(100),
        border = 'white',
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10))
#top_annotation = top_anno)


#matrix_1 <- RenameIdents(matrix,
#                           "0"="Tem",
#                           "1"="Tem",
#                           "2"="Naive CD8T",
#                           "3"="Tex",
#                           "4"="Naive CD8T",
#                           "5"="Tex",
#                           "6"="Tem",
#                           "7"="DEFB1+CD8T","8"="prolifer T","9"=""
#)
matrix_1 <- RenameIdents(matrix,"2"="CD8T","3"="CD8T","4"="CD8T","5"="CD8T","15"="CD8T")
matrix_1$celltype_sub <-Idents(matrix_1) 
DimPlot(matrix_1,reduction = "umap",label = T,group.by = "celltype_sub",split.by = "group")
pdf("macro_markers_dotplot.pdf",height = 3,width = 8)
DotPlot(matrix_1,features = c("GZMK","CXCR4","CXCR3","CD44","CCR7","SELL","PDCD1","LAG3","CTLA4","TIGIT","DEFB1","UGT2B7","MKI67","TOP2A"),assay = "RNA")
dev.off()
macro_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top5 <- macro_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf("matrix_top5heatmap.pdf",height = 10,width = 10)
DoHeatmap(matrix_1, features = top5$gene) + NoLegend()
dev.off()
pdf("annotation_1st.pdf",height = 10,width = 15)
DimPlot(matrix_1,reduction = "umap",label = T)
dev.off()

library(RColorBrewer)
colorcount <- length(unique(RCC_merge_filter@meta.data$celltype))
getPalette <- colorRampPalette(brewer.pal(12,"Paired"))
celltype_colors <- getPalette(colorcount)
DimPlot(RCC_merge_filter, reduction = "tsne",label = T,cols = celltype_colors)

CD8Tncc.markers <- FindMarkers(matrix, ident.1 = "nccRCC",group.by = "group", min.pct = 0.25,min.diff.pct = 0.25)

Cellratio <- prop.table(table(Idents(matrix_1), matrix_1$group), margin = 2)
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
pdf("macro_cellratio_cc_vs_ncc.pdf")
ggplot(Cellratio,aes(x =Var2, y= Freq, fill = Var1)) + 
  geom_bar(stat = "identity",position = "dodge")+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
sampleratio <- prop.table(table(matrix_1$orig.ident,Idents(matrix_1)), margin = 2)
sampleratio
sampleratio <- as.data.frame(sampleratio)
colourCount = length(unique(sampleratio$Var1))
library(ggplot2)
pdf("macro_sampleratio.pdf")
ggplot(sampleratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
Cellratio <- prop.table(table(Idents(matrix_1), matrix_1$orig.ident), margin = 2)
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
pdf("cellratio.pdf")
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  #coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

data <- GetAssayData(matrix_1, assay = 'RNA', slot = 'counts')
cell_metadata <- matrix_1@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


cds <- preprocess_cds(cds, num_dim = 100)


cds <- reduce_dimension(cds)

plot_cells(cds)

plot_cells(cds, color_cells_by="celltype_sub")
plot_cells(cds, color_cells_by="orig.ident")

cds <- align_cds(cds, alignment_group = "orig.ident")

cds <- reduce_dimension(cds)
#plot_cells(cds, genes=c("LAG3"))

#?reduce_dimension
#cds <- reduce_dimension(cds,umap.fast_sgd=T,cores = 3,reduction_method="tSNE")

#plot_cells(cds, reduction_method="tSNE")

#plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype_sub")


#plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE,reduction_method="tSNE")

#cds = align_cds(cds, num_dim = 100, alignment_group ="orig.ident")
#cds = reduce_dimension(cds,reduction_method="tSNE")
#plot_cells(cds, color_cells_by="orig.ident", label_cell_groups=FALSE,,reduction_method="tSNE")

cds = cluster_cells(cds, resolution=1e-5)
plot_cells(cds)

table(cds@clusters$UMAP$clusters)
cds@clusters$UMAP$partitions
#plot_cells(cds, color_cells_by="partition", group_cells_by="partition")

#plot_cells(cds, color_cells_by="celltype")

#marker_test_res <- top_markers(cds, group_cells_by="partition", 
reference_cells=1000, cores=8)


#top_specific_markers <- marker_test_res %>%
filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

#top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="celltype_sub",
                    ordering_type="maximal_on_diag",
                    max.size=3)





#top_specific_markers = marker_test_res %>%
#  filter(fraction_expressing >= 0.10) %>%
#  group_by(cell_group) %>%
#  top_n(3, pseudo_R2)

#top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))

#plot_genes_by_group(cds,
#                    top_specific_marker_ids,
#                    group_cells_by="celltype",
#                    ordering_type="cluster_row_col",
#                    max.size=3)

# 先将partitions的分组由因子型转为字符型
#colData(cds)$assigned_cell_type <- as.character(partitions(cds))

cds_subset <- choose_cells(cds)


cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "celltype_sub",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size=4,cell_size=1.5)


cds = order_cells(cds)

plot_cells(cds,
           color_cells_by = "group",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)
save_monocle_objects(cds=cds, directory_path='CD8T_1_cds_objects', comment='This is my macro cds. Stored 2023-04-09.')
cds <- load_monocle_objects(directory_path='CD8T_1_cds_objects')


table(matrix@meta.data$celltype)
# a helper function to identify the root principal points:
get_earliest_principal_node  <- function(cds, time_bin="T_cells"){
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

plot_cells(cds,
           color_cells_by = "group",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)


cds_sub <- choose_graph_segments(cds)




#Working with 3D trajectories

cds_3d = reduce_dimension(cds, max_components = 3)
cds_3d = cluster_cells(cds_3d)
cds_3d = learn_graph(cds_3d)
cds_3d = order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj = plot_cells_3d(cds_3d, color_cells_by="celltype")
cds_3d_plot_obj


ciliated_genes = top_specific_markers$gene_id[5:10]
cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

gene_fits = fit_models(cds_subset, model_formula_str = "~seurat_clusters")
fit_coefs
fit_coefs = coefficient_table(gene_fits)
# 挑出时间相关的组分
emb_time_terms = fit_coefs %>% filter(term == "seurat_clusters1")
emb_time_terms

emb_time_terms = fit_coefs %>% filter(term == "seurat_clusters1")

emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

plot_genes_violin(cds_subset[,], group_cells_by="celltype", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


neurons_cds <- cds[,grepl("Macrophage", colData(cds)$celltype, ignore.case=TRUE)]
plot_cells(neurons_cds, color_cells_by="partition")

pr_graph_test_res <- graph_test(neurons_cds, neighbor_graph="knn", cores=3)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df = find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)

cell_group_df = tibble::tibble(cell=row.names(colData(neurons_cds)), cell_group=neurons_cds@colData@listData[["seurat_clusters"]])
agg_mat = aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) = stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)


plot_cells(neurons_cds,
           genes=gene_module_df %>% filter(module %in% c(16,38,33,42)),
           group_cells_by="cluster",
           color_cells_by="cluster",
           show_trajectory_graph=FALSE)


plot_cells(cds,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


##找到影响发育轨迹的基因
trace('calculateLW', edit = T, where = asNamespace("monocle3"))
## change Matrix::rBind to rbind on line 93
# 使用neighbor_graph="principal_graph"来检验轨迹相邻的细胞的表达是否相关
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)



pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))

AFD_genes <-pr_deg_ids[1:10]
AFD_lineage_cds <- cds[rowData(cds)$gene_short_name %in% AFD_genes,colData(cds)$celltype_sub %in% c("Tex")
]
plot_genes_in_pseudotime(cds["PDCD1",],
                         color_cells_by = "celltype_sub",
                         min_expr=0.5)


cds_subset <- choose_cells(cds)

subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)

agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

plot_



dependencies ‘biomaRt’, ‘intrinsicDimension’, ‘parallelDist’, ‘R.cache’, ‘Signac’, ‘slingshot’ are not available for package ‘SCP’