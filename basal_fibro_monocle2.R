pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
library(pacman)
library(monocle)
basal_sub <- subset(epi_sub_3,celltype_sub=="basal_Krt15high"|celltype_sub=="basal_Krt15low")
fibro_sub <- subset(mesenchymal_sub_2,celltype_sub=="Fibro_Gpx3"|celltype_sub=="Fibro_Dpp4"|celltype_sub=="Fibro_Myoc")
basal_sub <- subset(basal_sub,downsample=500)
fibro_sub <- subset(fibro_sub,downsample=500)

basal_fibro_merge <- merge(basal_sub,fibro_sub)
basal_fibro_merge <- subset(mac_sub,downsample=500)
expr_matrix <- as(as.matrix(basal_fibro_merge@assays$RNA@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <-basal_fibro_merge@meta.data 
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(basal_fibro_merge),row.names = row.names(basal_fibro_merge))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)
#构建basal_fibro_cds_1对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)


basal_fibro_cds_1 <- newCellDataSet(expr_matrix,
                                    phenoData = pd,
                                    featureData = fd,
                                    lowerDetectionLimit = 1,
                                    expressionFamily = negbinomial.size())
basal_fibro_cds_1 <- estimateSizeFactors(basal_fibro_cds_1)
basal_fibro_cds_1 <- estimateDispersions(basal_fibro_cds_1)

basal_fibro_cds_1 <- detectGenes(basal_fibro_cds_1, min_expr = 0.125 ) #min_expr=0.125
print(head(fData(basal_fibro_cds_1)))
expressed_genes <- row.names(subset(fData(basal_fibro_cds_1), num_cells_expressed >= 10))
#expressed_genes <- row.names(basal_fibro_cds_1)
print(head(pData(basal_fibro_cds_1)))
pie <- ggplot(pData(basal_fibro_cds_1),
              aes(x = factor(1), fill = factor(seurat_clusters))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
diff_test_res <- differentialGeneTest(basal_fibro_cds_1[expressed_genes,],
                                      fullModelFormulaStr = "~celltype_sub")#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
ordering_genes <- subset(diff_test_res, qval < 0.01) 
ordering_genes <- ordering_genes[order(ordering_genes$qval,decreasing=F),]
ordering_genes <- rownames(ordering_genes)
basal_fibro_cds_1 <- setOrderingFilter(basal_fibro_cds_1, ordering_genes)
plot_ordering_genes(basal_fibro_cds_1)#出的图黑色的点表示用来构建轨迹的差异基因，
#灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)
basal_fibro_cds_1 <- reduceDimension(basal_fibro_cds_1, max_components = 2,
                                     method = 'DDRTree')
if(F){
  root_state = NULL
  num_paths = NULL
  reverse = NULL
  root_cell <- select_root_cell(basal_fibro_cds_1, root_state, reverse)
  basal_fibro_cds_1@auxOrderingData <- new.env(hash = TRUE)
  
  if (basal_fibro_cds_1@dim_reduce_type == "DDRTree") {
    if (is.null(num_paths) == FALSE) {
      message("Warning: num_paths only valid for method 'ICA' in reduceDimension()")
    }
    cc_ordering <- extract_ddrtree_ordering(basal_fibro_cds_1, root_cell)
    pData(basal_fibro_cds_1)$Pseudotime <- cc_ordering[row.names(pData(basal_fibro_cds_1)), ]$pseudo_time
    K_old <- reducedDimK(basal_fibro_cds_1)
    old_dp <- cellPairwiseDistances(basal_fibro_cds_1)
    old_mst <- minSpanningTree(basal_fibro_cds_1)
    old_A <- reducedDimA(basal_fibro_cds_1)
    old_W <- reducedDimW(basal_fibro_cds_1)
    basal_fibro_cds_1 <- project2MST(basal_fibro_cds_1, project_point_to_line_segment)
    minSpanningTree(basal_fibro_cds_1) <- basal_fibro_cds_1@auxOrderingData[[basal_fibro_cds_1@dim_reduce_type]]$pr_graph_cell_proj_tree
    root_cell_idx <- which(V(old_mst)$name == root_cell, arr.ind = T)
    cells_mapped_to_graph_root <- which(basal_fibro_cds_1@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex == root_cell_idx)
    if (length(cells_mapped_to_graph_root) == 0) {
      cells_mapped_to_graph_root <- root_cell_idx
    }
    cells_mapped_to_graph_root <- V(minSpanningTree(basal_fibro_cds_1))[cells_mapped_to_graph_root]$name
    tip_leaves <- names(which(degree(minSpanningTree(basal_fibro_cds_1)) == 1))
    root_cell <- cells_mapped_to_graph_root[cells_mapped_to_graph_root %in% tip_leaves][1]
    if (is.na(root_cell)) {
      root_cell <- select_root_cell(basal_fibro_cds_1, root_state, reverse)
    }
    basal_fibro_cds_1@auxOrderingData[[basal_fibro_cds_1@dim_reduce_type]]$root_cell <- root_cell
    cc_ordering_new_pseudotime <- extract_ddrtree_ordering(basal_fibro_cds_1, root_cell)
    pData(basal_fibro_cds_1)$Pseudotime <- cc_ordering_new_pseudotime[row.names(pData(basal_fibro_cds_1)), ]$pseudo_time
    if (is.null(root_state) == TRUE) {
      closest_vertex <- basal_fibro_cds_1@auxOrderingData[["DDRTree"]]$pr_graph_cell_proj_closest_vertex
      pData(basal_fibro_cds_1)$State <- cc_ordering[closest_vertex[, 1], ]$cell_state
    }
  }
}
basal_fibro_cds_1 <- orderCells(basal_fibro_cds_1)
save(basal_fibro_cds_1,file = "/home/liyang/new/basal_fibro_cds.rda")
plot_cell_trajectory(basal_fibro_cds_1, color_by = "celltype_sub",show_branch_points = F)+scale_color_manual(values =c("#FDB462", "#B3DE69","#A6CEE3" ,"#1F78B4","#33A02C"))+theme(legend.position = "right")
plot_cell_trajectory(basal_fibro_cds_1, color_by = "Pseudotime")+scale_color_viridis_c()+theme(legend.position = "right")#+scale_color_manual(values = brewer.pal("plasma"))
plot_cell_trajectory(basal_fibro_cds_1, color_by = "group")+facet_wrap("~group")+scale_color_manual(values = c("#FB9A99","#A6CEE3"))
plot_cell_trajectory(basal_fibro_cds_1, color_by = "group")# + facet_wrap("~orig.ident", nrow = 1)+scale_color_manual(values = brewer.pal(5,"Paired"))#以细胞状态上色（拆分）“分面”轨迹图，以便更容易地查看每个状态的位置。
plot_complex_cell_trajectory(basal_fibro_cds_1, x = 1, y = 2,
                             color_by = "seurat_clusters")+
  scale_color_manual(values = colour) +
  theme(legend.title = element_blank())#tree photo
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080")


sig_gene_names <- c("Tgm4","Tcaf3")
plot_pseudotime_heatmap(basal_fibro_cds_1[sig_gene_names,],
                        num_clusters = 1,
                        cores = 1,
                        show_rownames = T)

##选择前4个top基因并将其对象取出
keygenes <- head(ordering_genes,4)
basal_fibro_cds_1_subset <- basal_fibro_cds_1[keygenes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(basal_fibro_cds_1_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(basal_fibro_cds_1_subset, color_by = "celltype")
monocle::plot_genes_in_pseudotime(basal_fibro_cds_1[c("Cdh1","Cdh2","Vim","Fn1"),], color_by = "group")+scale_color_manual(values = c("#FB9A99","#A6CEE3"))
#指定基因
s.genes <- c("SELL","CCR7","IL7R", "CD84","CCL5","S100A4")
p1 <- plot_genes_jitter(basal_fibro_cds_1[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(basal_fibro_cds_1[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(basal_fibro_cds_1[s.genes,], color_by = "State")
monocle::plot_genes_in_pseudotime(basal_fibro_cds_1[c("Cdh1","Cdh2","Fn1","Vim"),], color_by = "orig.ident")
c("Ccl8","Ly6e","Ccl12","Maf")
c("Mndal","Iigp1","Pde7b","Ifi213")
c("Mmp12","Cd36")
#寻找拟时差异基因（qvalue体现基因与拟时的密切程度）绘制热图
#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(basal_fibro_cds_1[ordering_genes,], cores = 8,fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()
plot_pseudotime_heatmap(basal_fibro_cds_1[Time_genes,], num_clusters=3, show_rownames=F, return_heatmap=T)
#top100gene
Time_genes <- top_n(Time_diff, n = 100, desc(qval)) %>% pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(basal_fibro_cds_1[Time_genes,], num_clusters=3, show_rownames=F, return_heatmap=T)
save(basal_fibro_cds_1,Time_diff,file=".rda")
BEAM_res <- BEAM(basal_fibro_cds_1[ordering_genes,], branch_point = 1, cores = 8, progenitor_method = "duplicate") 
#这里用的是ordergene，也就是第六步dpFeature找出来的基因。如果前面用的是seurat的marker基因，记得改成express_genes
BEAM_res_1 <- BEAM(basal_fibro_cds_1, branch_point = 1, progenitor_method = "duplicate",cores = 8) #对2829个基因进行排序，运行慢
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
write.csv(BEAM_res, "/home/liyang/new/BEAM_res.csv", row.names = F)
tmp2 <- monocle::plot_genes_branched_heatmap(basal_fibro_cds_1[BEAM_res$gene_short_name[BEAM_res$qval < 1e-4],],
                                             branch_point = 1, #绘制的是哪个分支
                                             num_clusters = 3, #分
                                             use_gene_short_name = F,
                                             show_rownames = F,return_heatmap = T)


BEAM_genes <- subset(BEAM_res,
                     qval < 1e-4) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(basal_fibro_cds_1[BEAM_genes,], 
                                 num_clusters = 3, show_rownames = F, return_heatmap = T)
branchgene <- tmp2$annotation_row
