pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(pacman)
library(monocle)
expr_matrix <- as(as.matrix(luminal_sub_annotation@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <-luminal_sub_annotation@meta.data 
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(luminal_sub_annotation),row.names = row.names(luminal_sub_annotation))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)
#构建luminal_monocle对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)


luminal_monocle <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 1,
                      expressionFamily = negbinomial.size())
luminal_monocle <- estimateSizeFactors(luminal_monocle)
luminal_monocle <- estimateDispersions(luminal_monocle)

luminal_monocle <- detectGenes(luminal_monocle, min_expr = 0.125 ) #min_expr=0.125
print(head(fData(luminal_monocle)))
expressed_genes <- row.names(subset(fData(luminal_monocle), num_cells_expressed >= 10))
expressed_genes <- row.names(luminal_monocle)
print(head(pData(luminal_monocle)))
pie <- ggplot(pData(luminal_monocle),
              aes(x = factor(1), fill = factor(seurat_clusters))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
diff_test_res <- differentialGeneTest(luminal_monocle[expressed_genes,],
                                      fullModelFormulaStr = "~percent.mt",cores = 110)#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
ordering_genes <- subset(diff_test_res, qval < 0.01) 
ordering_genes <- ordering_genes[order(ordering_genes$qval,decreasing=F),]
ordering_genes <- rownames(ordering_genes)
luminal_monocle <- setOrderingFilter(luminal_monocle, ordering_genes)
plot_ordering_genes(luminal_monocle)#出的图黑色的点表示用来构建轨迹的差异基因，
#灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)
luminal_monocle <- reduceDimension(luminal_monocle, max_components = 2,
                       method = 'DDRTree')
if(F){
  root_state = NULL
  num_paths = NULL
  reverse = NULL
  root_cell <- select_root_cell(luminal_monocle, root_state, reverse)
  luminal_monocle@auxOrderingData <- new.env(hash = TRUE)
  
  if (luminal_monocle@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(luminal_monocle, root_cell)
    pData(luminal_monocle)$Pseudotime <- cc_ordering[row.names(pData(luminal_monocle)), ]$pseudo_time
    K_old <- reducedDimK(luminal_monocle)
    old_dp <- cellPairwiseDistances(luminal_monocle)
    old_mst <- minSpanningTree(luminal_monocle)
    old_A <- reducedDimA(luminal_monocle)
    old_W <- reducedDimW(luminal_monocle)
    luminal_monocle <- project2MST(luminal_monocle, project_point_to_line_segment)
    minSpanningTree(luminal_monocle) <- luminal_monocle@auxOrderingData[[luminal_monocle@dim_reduce_type]]$pr_graph_cell_proj_tree
    root_cell_idx <- which(V(old_mst)$name == root_cell, arr.ind = T)
    cells_mapped_to_graph_root <- which(luminal_monocle@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }
    cells_mapped_to_graph_root <- V(minSpanningTree(luminal_monocle))[cells_mapped_to_graph_root]$name
    tip_leaves <- names(which(degree(minSpanningTree(luminal_monocle)) == 1))
    root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
    if (is.na(root_cell)) {
      root_cell <- select_root_cell(luminal_monocle, root_state, reverse)
    }
    luminal_monocle@auxOrderingData[[luminal_monocle@dim_reduce_type]]$root_cell <- root_cell
    cc_ordering_new_pseudotime <- extract_ddrtree_ordering(luminal_monocle, root_cell)
    pData(luminal_monocle)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(luminal_monocle)), ]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- luminal_monocle@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      pData(luminal_monocle)$State <- cc_ordering[closest_vertex[, 1], ]$cell_state
    }
  }
}
luminal_monocle <- orderCells(luminal_monocle)
save(luminal_monocle,file = "luminal_monocle.rda")
plot_cell_trajectory(luminal_monocle, color_by = "seurat_clusters")+scale_color_manual(values = brewer.pal(5,"Paired"))
plot_cell_trajectory(luminal_monocle, color_by = "Pseudotime")
plot_cell_trajectory(luminal_monocle, color_by = "State")+scale_color_manual(values = c("#98FB98","#20B2AA","#FFA500"))
plot_cell_trajectory(luminal_monocle, color_by = "celltype_sub") + facet_wrap("~orig.ident", nrow = 1)+scale_color_manual(values = brewer.pal(5,"Paired"))#以细胞状态上色（拆分）“分面”轨迹图，以便更容易地查看每个状态的位置。
plot_complex_cell_trajectory(luminal_monocle, x = 1, y = 2,
                             color_by = "seurat_clusters")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank())#tree photo
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")


sig_gene_names <- c("Malat1","Nkx3-1")
plot_pseudotime_heatmap(luminal_monocle[sig_gene_names,],
                        num_clusters = 1,
                        cores = 1,
                        show_rownames = T)

##选择前4个top基因并将其对象取出
keygenes <- head(ordering_genes,4)
luminal_monocle_subset <- luminal_monocle[keygenes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(luminal_monocle_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(luminal_monocle_subset, color_by = "celltype")
plot_genes_in_pseudotime(luminal_monocle["Tgm4",], color_by = "celltype_sub")
#指定基因
s.genes <- c("SELL","CCR7","IL7R", "CD84","CCL5","S100A4")
p1 <- plot_genes_jitter(luminal_monocle[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(luminal_monocle[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(luminal_monocle[s.genes,], color_by = "State")
monocle::plot_genes_in_pseudotime(luminal_monocle[c("Malat1","Nkx3-1"),], color_by = "orig.ident")
c("Ccl8","Ly6e","Ccl12","Maf")
c("Mndal","Iigp1","Pde7b","Ifi213")
c("Mmp12","Cd36")
#寻找拟时差异基因（qvalue体现基因与拟时的密切程度）绘制热图
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(luminal_monocle[ordering_genes,], cores = 110,fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
p <- monocle::plot_pseudotime_heatmap(luminal_monocle[Time_genes,], num_clusters=3, show_rownames=F, return_heatmap=T)
#top100gene
Time_genes <- top_n(Time_diff, n = 2000, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(luminal_monocle[Time_genes,], cores = 110,num_clusters=3, show_rownames=F, return_heatmap=T)
clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
save(luminal_monocle,Time_diff,file=".rda")
BEAM_res <- BEAM(luminal_monocle[ordering_genes,], branch_point = 1, cores = 110, progenitor_method = "duplicate") 
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的marker基因，记得改成express_genes
BEAM_res <- BEAM(luminal_monocle, branch_point = 1, progenitor_method = "duplicate",cores = 110) #对2829个基因进行排序，运行慢
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
head(BEAM_res)
#           gene_short_name         pval         qval
# CD79A               CD79A 2.757676e-73 7.782161e-70
# TCL1A               TCL1A 1.574889e-65 2.222168e-62
# IGLL5               IGLL5 2.356778e-64 2.216942e-61
# S100A9             S100A9 1.504319e-58 1.061297e-55
# S100A8             S100A8 6.028175e-57 3.402302e-54
# LINC00926       LINC00926 3.180527e-55 1.495908e-52
write.csv(BEAM_res, "BEAM_res.csv", row.names = F)
tmp2 <- plot_genes_branched_heatmap(luminal_monocle[row.names(subset(BEAM_res,
                                                 qval < 1e-4)),],
                            branch_point = 1, #绘制的是哪个分支
                            num_clusters = 4, #分成几个cluster，根据需要调整
                            cores = 1,
                            use_gene_short_name = F,
                            show_rownames = F,return_heatmap = T)


BEAM_genes <- subset(BEAM_res,
                     qval < 1e-4) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(luminal_monocle[BEAM_genes,],  branch_point = 1,
                                 num_clusters = 3, show_rownames = F, return_heatmap = T)
branchgene <- tmp2$annotation_row
