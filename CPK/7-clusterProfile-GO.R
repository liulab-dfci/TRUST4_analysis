## cluster profile for GO analysis
setwd("/liulab/zouyang/TCGA/IGH/6_trust4_clustering/TPM_corr_CPK/")
#BiocManager::install("clusterProfiler")
options(digits=2) # for the number of digits to show
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)# for shorten the long-word GO name
# dotplot() + scale_y_discrete(labels=function(x) str_wrap(x, width=10))
my_theme = theme( legend.title = element_text(size = 15),
                 legend.text = element_text(size = 15),
                 legend.key.width=unit(1, "lines"),
                 plot.title = element_text(size = 15, face="bold", hjust=0.5),
                 #plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),
                 axis.text.x = element_text(angle = 0, size = 15,face="bold",
                                            color = "black"),
                 axis.title.x=element_text(face = "bold"),
                 axis.text.y = element_text(face="bold", size = 15),
                 #panel.background = element_rect(fill = NA, colour = "black"),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 panel.grid.major = element_line(linetype="dashed"),
                 panel.grid.minor = element_line(linetype="dashed"))#
#==============Ig fraction================
# IGH go
IGHgene = read.table("top200_genes_in_CPKcomplete.IGH.txt")
IGHego <- enrichGO(gene          = IGHgene$V1,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(IGHego))
IGHego.sim <- simplify( IGHego, cutoff = 0.7, by = "p.adjust", select_fun = min)
pdf("Fig_GO_of_correlated_genes_with_CPK.pdf",width = 9)
#dotplot(IGHego,showCategory=10)
dotplot(IGHego.sim,showCategory=10) + scale_y_discrete(labels=function(x) str_wrap(x, width=45)) +
  theme_classic2() + my_theme

# TRB GO
TRBgene = read.table("top200_genes_in_CPKcomplete.TRB.txt")
TRBego <- enrichGO(gene       = TRBgene$V1,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)
head(summary(TRBego))
TRBego.sim <- simplify( TRBego, cutoff = 0.7, by = "p.adjust", select_fun = min)
#dotplot(TRBego,showCategory=10)
dotplot(TRBego.sim,showCategory=10) + scale_y_discrete(labels=function(x) str_wrap(x, width=45)) +
  theme_classic2() + my_theme
dev.off()


#============== (IgG1+IgG3)% ================
# IGH go
IgG13gene = read.table("top200_genes_in_IgG13_fraction.absolute.txt")
IgG13ego <- enrichGO(gene          = IgG13gene$V1,
                   OrgDb         = org.Hs.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
head( summary(IgG13ego) )
IgG13ego.sim <- simplify( IgG13ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

pdf("Fig_GO_of_correlated_genes_with_IgG13fraction.pdf",width = 9)
dotplot(IgG13ego,showCategory=10)
dotplot(IgG13ego.sim,showCategory=10) + scale_y_discrete(labels=function(x) str_wrap(x, width=45)) +
  theme_classic2() + my_theme
dev.off()

save( IGHego, IGHego.sim, TRBego, TRBego.sim, IgG13ego, IgG13ego.sim, file = "Simplify_GO_for_IGH_TRB_and_IgG13fraction.Rdata")
