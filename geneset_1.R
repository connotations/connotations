ROS <- read.csv("/home/liyang/liyang/prostate-age/GOBP_CELLULAR_RESPONSE_TO_REACTIVE_OXYGEN_SPECIES.v2023.1.Mm.tsv")
ROS <- ROS[168:313,]
ROS[1] <- "Adprs"
ROS <- intersect(ROS,rownames(epi_sub_3@assays$RNA))
ROS<- as.data.frame(ROS)
ROS <- as.list(ROS)
# 计算基因集得分
prostate_merge_annotation_6st <- AddModuleScore(object = prostate_merge_annotation_6st, 
                                                features =ROS,
                                                name = "ROS_score",
                                                ctrl = 100)
VlnPlot(object = prostate_merge_annotation_4st, 
        features = "ROS_score1",
        cols = c("pink",  "#048BA8"),group.by = "group",pt.size = 0,
        combine = FALSE)
# 绘制基因集得分的密度图
VlnPlot(object = prostate_merge_annotation_4st, 
        features = "telomere_score1",
        cols = c("pink",  "#048BA8"),group.by = "celltype",pt.size = 0,split.by = "group",
        combine = TRUE)
RidgePlot(prostate_merge_annotation_4st, features = "infla_resp_score1",group.by = 'group',cols = c("pink",  "#048BA8"))
myeloid_sub <- prostate
saveRDS(myeloid_sub,file = "prostate_merge_uptodate.rds")
senescence <- intersect(senescence,rownames(prostate_merge_annotation_6st@assays$RNA))
senescence <- as.data.frame(senescence )
senescence <- as.list(senescence )
prostate_merge_annotation_6st <- AddModuleScore(object =prostate_merge_annotation_6st, 
                                                features = senescence,
                                                name = "senescence_score",
                                                ctrl = 100)
senescence <- read.csv("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP.v2023.1.Mm.tsv")
senescence <- senescence[58:97,]
senescence <- intersect(senescence,rownames(prostate_merge_annotation_4st@assays$RNA))
senescence <- as.data.frame(senescence )
senescence <- as.list(senescence )
infla_resp <- read.csv("../HALLMARK_INFLAMMATORY_RESPONSE.v2023.1.Mm.tsv")
infla_resp <- infla_resp[226:422,]
infla_resp[1] <- "Mmp14"
infla_resp <- intersect(infla_resp,rownames(prostate_merge_annotation_4st@assays$RNA))
infla_resp<- as.data.frame(infla_resp )
infla_resp <- as.list(infla_resp )
T_sub_2<- AddModuleScore(object =T_sub_2, 
                         features = infla_resp,
                         name = "infla_resp_score",
                         ctrl = 100)
chronic_infla_resp <- read.csv("/home/liyang/liyang/prostate-age/GOBP_CHRONIC_INFLAMMATORY_RESPONSE.v2023.1.Mm.tsv")
chronic_infla_resp <- chronic_infla_resp[46:67,]
chronic_infla_resp[1] <- "Adora2b"
chronic_infla_resp <- intersect(chronic_infla_resp,rownames(epi_sub_3@assays$RNA))
chronic_infla_resp<- as.data.frame(chronic_infla_resp )
chronic_infla_resp <- as.list(chronic_infla_resp )
epi_sub_3<- AddModuleScore(object =epi_sub_3, 
                           features = chronic_infla_resp,
                           name = "chronic_infla_resp_score",
                           ctrl = 100)

telomere <- read.csv("../REACTOME_TELOMERE_MAINTENANCE.v2023.1.Mm.tsv")
telomere <- telomere[65:113,]
telomere[1] <- "Rpa1"
telomere <- intersect(telomere,rownames(prostate_merge_annotation_4st@assays$RNA))
telomere<- as.data.frame(telomere )
telomere <- as.list(telomere )
prostate_merge_annotation_4st <- AddModuleScore(object =prostate_merge_annotation_4st, 
                                                features = telomere,
                                                name = "telomere_score",
                                                ctrl = 100)

DNArepair <- read.csv("HALLMARK_DNA_REPAIR.v2023.1.Mm.tsv")
DNArepair <- DNArepair[177:324,]
DNArepair[1] <- "Lig1"
DNArepair <- intersect(DNArepair,rownames(prostate_merge_annotation_4st@assays$RNA))
DNArepair<- as.data.frame(DNArepair)
DNArepair <- as.list(DNArepair)
prostate_merge_annotation_4st <- AddModuleScore(object =prostate_merge_annotation_4st, 
                                                features =DNArepair,
                                                name = "DNArepair_score",
                                                ctrl = 100)
regeneration <- read.csv("GOBP_REGENERATION.v2023.1.Mm.tsv")
View(regeneration)
regeneration <- regeneration[182:343,]
regeneration[1] <- "Lnpp5f"
regeneration <- intersect(regeneration,rownames(prostate_merge_annotation_4st@assays$RNA))
regeneration<- as.data.frame(regeneration)
regeneration <- as.list(regeneration)
prostate_merge_annotation_4st <- AddModuleScore(object =prostate_merge_annotation_4st, 
                                                features =regeneration,
                                                name = "regeneration_score",
                                                ctrl = 100)
#ROS
ROS_score <- FetchData(prostate_merge_annotation_6st,c("ROS_score1","group","celltype"))
myeloid_ROS_score <- FetchData(myeloid_sub_1,c("group","ROS_score1"))
B_ROS_score <- FetchData(B_sub,c("group","ROS_score1"))
T_ROS_score <- FetchData(T_sub,c("group","ROS_score1"))
ggboxplot(ROS_score, x = "celltype", y = "ROS_score1",
          # 配色方案 ?ggboxplot
          color = "group", palette = c("#FB9A99","#A6CEE3" )
)#+ stat_compare_means(label =  "p.signif", label.x = 1.5)
#  Add p-value
p + stat_compare_means(label =  "p.signif", label.x = 1.5) #default Wilcoxon
p + stat_compare_means(method = "t.test")
compare_means(ROS_score1~group,data = ROS_score)
ggplot(ROS_score,aes(x=group,y=ROS_score1,color=group))+facet_wrap("~celltype")+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)+scale_color_manual(values = c("#FB9A99","#A6CEE3" ))
ggplot(T_ROS_score,aes(x=group,y=ROS_score1,color=group))+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(ROS_score,aes(x=ROS_score1,color=group,fill=group))+ geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1) +scale_fill_manual(values = c("#FB9A99","#A6CEE3" ))+scale_color_manual(values = c("#FB9A99","#A6CEE3" ))+
  theme_classic()+geom_vline(aes(xintercept=median(ROS_score$ROS_score1[ROS_score$group=="Aged"])),linetype=5,col="#FB9A99")+geom_vline(aes(xintercept=median(ROS_score$ROS_score1[ROS_score$group=="Young"])),linetype=5,col="#A6CEE3")
#infla resp
chronic_infla_resp_score <- FetchData(prostate_merge_annotation_6st,c("group","chronic_infla_resp_score1","celltype_sub"))
B_infla_resp_score <- FetchData(B_sub,c("group","infla_resp_score1"))
T_infla_resp_score <- FetchData(T_sub_2,c("group","infla_resp_score1","celltype_sub"))
ggboxplot(chronic_infla_resp_score, x = "group", y = "chronic_infla_resp_score1",
          # 配色方案 ?ggboxplot
          color = "group", palette = c("#FB9A99","#A6CEE3" )
)+ stat_compare_means(label =  "p.signif", label.x = 1.5)
#  Add p-value
p + stat_compare_means(label =  "p.signif", label.x = 1.5) #default Wilcoxon
p + stat_compare_means(method = "t.test")
compare_means(ROS_score1~group,data = ROS_score)
ggplot(chronic_infla_resp_score,aes(x=group,y=chronic_infla_resp_score1,color=group))+facet_wrap("~celltype_sub")+geom_boxplot(width=0.1,cex=1.2)+theme_cowplot()+ stat_compare_means(label =  "p.signif", label.x = 1.5)+theme_classic()
ggplot(T_infla_resp_score,aes(x=group,y=infla_resp_score1,color=group))+facet_wrap ("~celltype_sub")+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)

ggplot(chronic_infla_resp_score,aes(x=chronic_infla_resp_score1,color=group,fill=group))+ geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1) +scale_fill_manual(values = c("#FB9A99","#A6CEE3" ))+scale_color_manual(values = c("#FB9A99","#A6CEE3" ))+
  theme_classic()+geom_vline(aes(xintercept=median(chronic_infla_resp_score$chronic_infla_resp_score1[chronic_infla_resp_score$group=="Aged"])),linetype=5,col="#FB9A99")+geom_vline(aes(xintercept=median(chronic_infla_resp_score$infla_resp_score1[chronic_infla_resp_score$group=="Young"])),linetype=5,col="#A6CEE3")
#senescence
senescence_score <- FetchData(prostate_merge_annotation_6st,c("senescence_score1","group","celltype"))
luminal_senescence_score <- FetchData(prostate_merge_annotation_4st[,prostate_merge_annotation_4st$celltype=="luminal"],c("group","senescence_score1"))
B_senescence_score <- FetchData(B_sub,c("group","senescence_score1"))
T_senescence_score <- FetchData(T_sub,c("group","senescence_score1"))
ggboxplot(senescence_score, x = "celltype", y = "senescence_score1",
          # 配色方案 ?ggboxplot
          color = "group", palette = c("#FB9A99","#A6CEE3" ))
#  Add p-value
p + stat_compare_means(label =  "p.signif", label.x = 1.5) #default Wilcoxon
p + stat_compare_means(method = "t.test")
compare_means(ROS_score1~group,data = ROS_score)
ggplot(senescence_score,aes(x=group,y=senescence_score1,color=group))+facet_wrap("~celltype")+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)+scale_color_manual(values = c("#FB9A99","#A6CEE3" ))
ggplot(B_senescence_score,aes(x=group,y=senescence_score1,color=group))+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)

ggplot(senescence_score,aes(x=senescence_score1,color=group,fill=group))+ geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1) +scale_fill_manual(values = c("#FB9A99","#A6CEE3" ))+scale_color_manual(values = c("#FB9A99","#A6CEE3" ))+
  theme_classic()+geom_vline(aes(xintercept=median(senescence_score$senescence_score1[senescence_score$group=="Aged"])),linetype=5,col="#FB9A99")+geom_vline(aes(xintercept=median(senescence_score$senescence_score1[senescence_score$group=="Young"])),linetype=5,col="#A6CEE3")
#DNArepair
DNArepair_score <- FetchData(prostate_merge_annotation_4st,c("DNArepair_score1","group"))
myeloid_DNArepair_score <- FetchData(myeloid_sub_2,c("group","DNArepair_score1"))
B_DNArepair_score <- FetchData(B_sub,c("group","DNArepair_score1"))
T_DNArepair_score <- FetchData(T_sub,c("group","DNArepair_score1"))
ggboxplot(ROS_score, x = "group", y = "ROS_score1",
          # 配色方案 ?ggboxplot
          color = "group", palette = "ROS_score1",
          add = "jitter",size = 0.5)

p + stat_compare_means(label =  "p.signif", label.x = 1.5) #default Wilcoxon
p + stat_compare_means(method = "t.test")
compare_means(ROS_score1~group,data = ROS_score)
ggplot(myeloid_DNArepair_score,aes(x=group,y=DNArepair_score1,color=group))+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(B_DNArepair_score,aes(x=group,y=DNArepair_score1,color=group))+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(DNArepair_score,aes(x=group,y=DNArepair_score1,color=group))+facet_wrap ("~celltype")+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_cowplot()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(DNArepair_score,aes(x=DNArepair_score1,color=group,fill=group))+ geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1) +scale_fill_manual(values = c("#FB9A99","#A6CEE3" ))+scale_color_manual(values = c("#FB9A99","#A6CEE3" ))+
  theme_classic()+geom_vline(aes(xintercept=median(DNArepair_score$DNArepair_score1[DNArepair_score$group=="Aged"])),linetype=5,col="#FB9A99")+geom_vline(aes(xintercept=median(DNArepair_score$DNArepair_score1[DNArepair_score$group=="Young"])),linetype=5,col="#A6CEE3")
#telomere
telomere_score <- FetchData(prostate_merge_annotation_4st,c("group","telomere_score1","celltype"))
myeloid_DNArepair_score <- FetchData(myeloid_sub_2,c("group","DNArepair_score1"))
B_DNArepair_score <- FetchData(B_sub,c("group","DNArepair_score1"))
T_DNArepair_score <- FetchData(T_sub,c("group","DNArepair_score1"))
ggboxplot(ROS_score, x = "group", y = "telomere_score1",
          # 配色方案 ?ggboxplot
          color = "group", palette = "ROS_score1",
          add = "jitter",size = 0.5)

p + stat_compare_means(label =  "p.signif", label.x = 1.5) #default Wilcoxon
p + stat_compare_means(method = "t.test")
compare_means(ROS_score1~group,data = ROS_score)
ggplot(telomere_score,aes(x=group,y=telomere_score1,color=group))+facet_wrap ("~celltype")+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(B_DNArepair_score,aes(x=group,y=DNArepair_score1,color=group))+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(DNArepair_score,aes(x=group,y=DNArepair_score1,color=group))+facet_wrap ("~celltype")+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_cowplot()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(telomere_score,aes(x=telomere_score1,color=group,fill=group))+ geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1) +
  theme_classic()+geom_vline(aes(xintercept=median(telomere_score$telomere_score1[telomere_score$group=="Aged"])),linetype=5,col="#FB9A99")+geom_vline(aes(xintercept=median(telomere_score$telomere_score1[telomere_score$group=="Young"])),linetype=5,col="#19D0D6")
#regeneration
regeneration_score <- FetchData(prostate_merge_annotation_4st,c("regeneration_score1","group","celltype"))
myeloid_DNArepair_score <- FetchData(myeloid_sub_2,c("group","DNArepair_score1"))
B_DNArepair_score <- FetchData(B_sub,c("group","DNArepair_score1"))
T_DNArepair_score <- FetchData(T_sub,c("group","DNArepair_score1"))
ggboxplot(ROS_score, x = "group", y = "telomere_score1",
          # 配色方案 ?ggboxplot
          color = "group", palette = "ROS_score1",
          add = "jitter",size = 0.5)

p + stat_compare_means(label =  "p.signif", label.x = 1.5) #default Wilcoxon
p + stat_compare_means(method = "t.test")
compare_means(ROS_score1~group,data = ROS_score)
ggplot(regeneration_score,aes(x=group,y=regeneration_score1,color=group))+facet_wrap ("~celltype")+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(B_DNArepair_score,aes(x=group,y=DNArepair_score1,color=group))+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_classic()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(regeneration_score,aes(x=group,y=regeneration_score1,color=group))+facet_wrap ("~celltype")+geom_violin()+geom_boxplot(width=0.1,cex=1.2)+theme_cowplot()+ stat_compare_means(label =  "p.signif", label.x = 1.5)
ggplot(regeneration_score,aes(x=regeneration_score1,color=group,fill=group))+ geom_density(adjust=1.5, alpha=0.5,lwd=1.5,linetype = 1) +scale_fill_manual(values = c("#FB9A99","#A6CEE3" ))+scale_color_manual(values = c("#FB9A99","#A6CEE3" ))+
  theme_classic()+geom_vline(aes(xintercept=median(regeneration_score$regeneration_score1[regeneration_score$group=="Aged"])),linetype=5,col="#FB9A99")+geom_vline(aes(xintercept=median(regeneration_score$regeneration_score1[regeneration_score$group=="Young"])),linetype=5,col="#A6CEE3")
