basal_Krt15low_degs_age <- FindMarkers(epi_sub_4[,epi_sub_4$celltype_sub=="basal_Krt15low"]
                                      ,ident.1 ="Aged"
                                      ,ident.2="Young",group.by = "group")

basal_Krt15low_degs_age$cluster <- "basal_Krt15low"
basal_Krt15low_degs_age$gene <- rownames(basal_Krt15low_degs_age)
