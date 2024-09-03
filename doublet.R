library(DoubletFinder)
sweep.res.list_2 <- paramSweep_v3(prostate_merge_2, PCs = 1:30, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list_2, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
DoubletRate = ncol(prostate_merge_2)*8*1e-6
homotypic.prop <- modelHomotypic(prostate_merge_2$celltype)
nExp_poi <- round(DoubletRate*ncol(prostate_merge_2))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
prostate_merge_2 <- doubletFinder_v3(prostate_merge_2, PCs = 1:30, pN = 0.25, pK = pK_bcmvn, 
                          nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
DimPlot(prostate_merge_2,reduction = "umap",group.by = "DF.classifications_0.25_0.001_1349",pt.size = 1)
prostate_merge_3 <- subset(prostate_merge_2,DF.classifications_0.25_0.001_1349=="Singlet")

