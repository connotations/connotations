Idents(TNK_ps) <- "sample"
OE_T_marker <- FindAllMarkers(TNK_ps)

OE_up <- OE_T_marker$gene[OE_T_marker$cluster=="OE"&OE_T_marker$avg_log2FC>0.5]
#lymphatic endo_go_down <- subset(lymphatic endo_gos,p_val_adj<0.05&avg_log2FC< -0.15)
OE_up <- enrichGO(gene          = OE_up,
                        #universe     = row.names(CD8Tncc._markers),
                        OrgDb         = 'org.Mm.eg.db',
                        keyType       = 'SYMBOL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
OE_up <- data.frame(OE_up)
View(OE_up)
OE_up$group <- "OE"
up <- OE_up[order(OE_up$p.adjust),]
up <- up[1 : 10,]

OE_down <- OE_T_marker$gene[OE_T_marker$cluster=="NC"&OE_T_marker$avg_log2FC>0.5]
#lymphatic endo_go_down <- subset(lymphatic endo_gos,p_val_adj<0.05&avg_log2FC< -0.15)
OE_down <- enrichGO(gene          = OE_down,
                  #universe     = row.names(CD8Tncc._markers),
                  OrgDb         = 'org.Mm.eg.db',
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05)
OE_down <- data.frame(OE_down)
View(OE_down)
OE_down$LogP <- log10(OE_down$p.adjust)
OE_down$group <- "NC"
down <- OE_down[order(OE_down$p.adjust),]
down <- down[1 : 10,]
down$Count <- -down$Count
A <- rbind(up,down)
  
library(ggplot2)
OE_up$LogP <- -log10(OE_up$p.adjust)
A <- OE_up
A <- A[order(A$p.adjust),]
A <- A[1 : 10,]
ggplot(A,aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=LogP)
           ,stat='identity')+
  scale_fill_gradient(low="#FFCC33",high="#CC6666")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
ggplot(A, aes(x = reorder(Description,Count), Count,fill=LogP)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6))+
  labs(x = "",
       y="Count")+
  scale_fill_gradient2(low="#CC6666",high="#8DD3C7")#设置颜色

intersect_down <- enrichGO(gene          = intersect(cluster2,rownames(plasma_age_markers[plasma_age_markers$p_val_adj<0.05&plasma_age_markers$avg_log2FC< -0.15,])),
                           #universe     = row.names(CD8Tncc._markers),
                           OrgDb         = 'org.Mm.eg.db',
                           keyType       = 'SYMBOL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.01,
                           qvalueCutoff  = 0.05)
intersect_down <- data.frame(intersect_down)
View(intersect_down)
intersect_down <- intersect_down[order(intersect_down$p.adjust),]
intersect_down_top10 <- intersect_down[1 : 8,]
ggplot(data=intersect_up, aes(x=Description,y=log10_p)) + 
  geom_bar(stat="identity", width=0.8,fill='#DC143C') + 
  coord_flip() +  xlab("GO terms") + ylab("-Log10P") + 
  theme_classic()
ggplot(data=intersect_down_top10, aes(x=Description,y=log10_p)) + 
  geom_bar(stat="identity", width=0.8,fill='#0000CD') + 
  coord_flip() +  xlab("GO terms") + ylab("-Log10P") + 
  theme_classic()+ scale_x_discrete(labels = function(x) str_wrap(x, width = 40) )
