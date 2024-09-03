library(ComplexHeatmap)
mat <- GetAssayData(Myeloid_filter,slot = "scale.data")
gene_features <- top20
cluster_info <- sort(Myeloid_filter$celltype_sub)
mat <- as.matrix(mat[top20$gene,names(cluster_info)])
gene <- c("Chil3","Ly6c2","Vcan","Cx3cr1","Jun","Cd74","Mki67","Top2a","Mmp12","Mmp13","Ifit3","Cmpk2","Ifit2","Ccl7","Ccl2","Ccl12","Lyve1","Lrg1")
gene_pos <- which(rownames(mat)%in%gene)
row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,labels = gene))

names(colour) <- levels(cluster_info)
top_anno <- HeatmapAnnotation(cluster=anno_block(gp=gpar(fill=colour),
                                                 labels = levels(cluster_info),
                                                 labels_gp = gpar(cex=0.5,col="white")))
library(circlize)
col_fun <- colorRamp2(c(-4,0,4),c("#377EB8","white","#E41A1C"))
Heatmap(mat,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,column_split = cluster_info,
        right_annotation = row_anno,
        top_annotation = top_anno,
        column_title = NULL,
        heatmap_legend_param = list(title="Expression",title_position="leftcenter-rot"),
        col = col_fun)
