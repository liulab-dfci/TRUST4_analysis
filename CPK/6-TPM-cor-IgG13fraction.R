library(dplyr)
library(data.table)
library(tidyr)
library(survcomp)
library(ggpubr)
library(ggplot2)
setwd('/liulab/zouyang/TCGA/IGH/6_trust4_clustering/TPM_corr_CPK/')
source("/liulab/zouyang/TCGA/IGH/NewCgene/SHMratio_res/heatmap.r")
my_theme = theme(legend.title = element_blank(),
                 legend.text = element_text(size=17, face = "bold"),
                 legend.key.width=unit(1, "lines"),
                 plot.title = element_text(size=10, face="bold", hjust=0.5),
                 #plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),
                 axis.text.x = element_text(size=17, face="bold"),
                 axis.title.x = element_text(size=14, face="bold"),
                 axis.text.y = element_text(size = 17, face="bold"),
                 axis.title.y = element_text(size = 14,face="bold"),
                 panel.grid.major = element_blank())

testMetric = "IgG13_fraction"
## correlation between CPK and gene expression
corr_path = './'
exprCPK = NULL
for(f in list.files(corr_path, pattern= 'corr_TPM_and_IgG13_fraction_spearman.*.Rdata')) {
    load(paste0(corr_path, f))
    if( length(exprPropIg.pcor)>1 ){
      cancerType = gsub('corr_TPM_and_IgG13_fraction_spearman.','',f)
      cancerType = gsub('.Rdata','',cancerType)
      exprCPK = rbind(exprCPK, data.frame(exprPropIg.pcor,Disease=cancerType,expr=rownames(exprPropIg.pcor)))
    }
}
exprCPK$padj = p.adjust(exprCPK$pvalue, method = "BH")
save(exprCPK,file="AllCorr_TPM_and_IgG13_fraction_spearman.Rdata")
exprCPKcor = spread(exprCPK[,c(1,3:4)],Disease,cor)
exprCPKpval = spread(exprCPK[,c(2:4)],Disease,pvalue)
exprCPKpval$combineP = apply(exprCPKpval[,-1],1,function(x) {combine.test(x,na.rm=T)})
exprCPKpval$padj = p.adjust(exprCPKpval$combineP, method='BH')
exprCPKpval = exprCPKpval[order(exprCPKpval$padj),]

## check Trust3 picked genes
t3bcell = c('CD37','FCGR3A','CD69')
t3oncogene = c('MYC','KRAS','MKI67')
t3tcell = c('CTLA4','PDCD1','CD40LG')
t3 = c(t3bcell,t3oncogene,t3tcell)
# CPK
which(exprCPKpval$expr %in% t3bcell)
which(exprCPKpval$expr %in% t3tcell)
which(exprCPKpval$expr %in% t3oncogene)

# top20 genes
cutoff=20
exprCPKcor$sign = apply(exprCPKcor[,-1],1,function(x) {length(which(sign(x)>0))})
exprCPKcor$medianCor = apply(exprCPKcor[,2:34],1,function(x) {median(x,na.rm=T)})
exprCPKpvalSort = merge(exprCPKpval,exprCPKcor[,c(1,35,36)],by='expr')
exprCPKpvalSort = exprCPKpvalSort[order(exprCPKpvalSort$padj,exprCPKpvalSort$medianCor),]
exprCPKpvalSort = exprCPKpvalSort[order(abs(exprCPKpvalSort$medianCor),decreasing=T),]
save(exprCPKpvalSort,file='geneOrder_IgG13_fraction.Rdata')
write.table(na.omit(exprCPKpvalSort[,c('expr','medianCor')]),"gene_corr_IgG13_fraction.rnk",sep='\t',row.names=F,quote=F)
write.table(na.omit(exprCPKpvalSort[,c('expr')]),"gene_corr_IgG13_fraction.gene.rnk",sep='\t',row.names=F,quote=F)

## rank gene by CPK
top200gene = exprCPKpvalSort[ order( abs(exprCPKpvalSort$medianCor), decreasing = T), ]
top200gene = as.character(top200gene$expr[1:200])
write.table( top200gene, "top200_genes_in_IgG13_fraction.absolute.txt", sep=",", quote=F,row.names=F,col.names=F)
exprCPKpvalSortSub = subset(exprCPKpvalSort, -log10(padj)>= 0 )
exprCPKpvalSortSub = exprCPKpvalSortSub[order(exprCPKpvalSortSub$padj,decreasing=F),]
exprCPKpvalSortSub = exprCPKpvalSortSub[order(exprCPKpvalSortSub$medianCor,decreasing=T),]
top10 = as.character(head(exprCPKpvalSortSub$expr,200))
IgG = subset(exprCPKpvalSortSub, medianCor>0.1)
#IgG.top20 = as.character(IgG$expr)
IgG.top20 = as.character(exprCPKpvalSortSub$expr[1:25])
exprCPKpvalSortSub = exprCPKpvalSortSub[order(exprCPKpvalSortSub$padj,decreasing=F),]
exprCPKpvalSortSub = exprCPKpvalSortSub[order(exprCPKpvalSortSub$medianCor),]
tail10 = as.character(head(exprCPKpvalSortSub$expr,200))
IgA.top20 = as.character(exprCPKpvalSortSub$expr[1:25])
IgA = subset(exprCPKpvalSortSub, medianCor< -0.1)
#IgA.top20 = as.character(IgA$expr)
write.table(IgG.top20,"top200_genes_in_IgG13_fraction.positive.txt",sep=",",quote=F,row.names=F,col.names=F)
write.table(IgA.top20,"top200_genes_in_IgG13_fraction.negative.txt",sep=",",quote=F,row.names=F,col.names=F)

exprCPKgroup <- exprCPK %>% group_by(expr) %>% summarise(corr=median(cor,na.rm=T),pvalue=combine.test(pvalue,na.rm=T))
exprCPKgroup$padj = p.adjust(exprCPKgroup$pvalue, method='BH')
exprCPKgroup = exprCPKgroup[order(exprCPKgroup$padj),]
exprCPKgroup$log10padj = -log10(exprCPKgroup$padj)
exprCPKgroup$labels = sapply(exprCPKgroup$expr, function(x) {
  if(x %in% IgG.top20){
    return('IgG-positive')
  }else if(x %in% IgA.top20){
    return('IgG-negative')
  }else{
    return('others')
  }
})

pdf('Fig_rank_gene_by_correlation_TPM_IgG13_fraction.pdf')
exprCPKgroup$labels <-factor(exprCPKgroup$labels,levels=c("IgG-positive","IgG-negative","others"))
ggscatter(exprCPKgroup,x='corr',y='log10padj',color='labels',palette=c(get_palette(palette='npg',2),"grey"),label='expr',repel=T,label.select=c(IgG.top20,IgA.top20),
      font.label = c(12, "bold"), size = 1.5,
  xlab="correlation between gene expression and IgA%",ylab="-log10(Padj)") + my_theme
  #geom_hline( yintercept = 50, linetype = "dashed") + my_theme
dev.off()

## heatmap of top 10 genes
exprCPK.sub = subset(exprCPK, expr %in% c(IgG.top20,IgA.top20))
exprCPK.sub$expr = factor(exprCPK.sub$expr, levels=c(IgG.top20,IgA.top20))
maxVal = round(max(exprCPK$cor,na.rm=T),2)
minVal = round(min(exprCPK$cor,na.rm=T),2)
exprCPK$cor[ exprCPK$cor > maxVal] = maxVal
exprCPK$cor[ exprCPK$cor < minVal] = minVal
p1<-corPlot(exprCPK.sub,'Disease','expr',rho='cor',"pvalue",max_value = maxVal,min_value = minVal, box_size = 5.5, bar_anno = c(minVal, 0 , maxVal))
ggsave(p1,file="Fig_heatmap_TPM_IgG13_fraction.pdf",width=8,height=8,dpi=1000)

##  label cytokine genes
cytokine = read.table("/homes/zouyang/Immunotherapy/src/cytokine_list.txt",header=F,stringsAsFactor=F)

exprCPKgroup$CtyLabels = sapply( as.character(exprCPKgroup$expr), function(x) {
  tmpdata = subset(exprCPKgroup, expr==x)
  if( x %in% cytokine$V1 ){
    if( !is.na(tmpdata$corr) & abs(tmpdata$corr)>0.1 ){
      return('cytokine')
    }else{
      #print( tmpdata$labels )
      return( tmpdata$labels)
    }
  }else{
    return( tmpdata$labels)
  }
})
table(exprCPKgroup$CtyLabels )
cytokinePick = subset(exprCPKgroup, CtyLabels != "others" )
cytokineCorr = subset(exprCPKgroup, gene %in% cytokine$V1)
summary(cytokineCorr$corr)
summary(cytokineCorr$padj)
pdf('Fig_rank_gene_by_correlation_TPM_IgG13_fraction.labelCytokine.pdf')
exprCPKgroup$CtyLabels <-factor(exprCPKgroup$CtyLabels,levels=c("IgG-positive","IgG-negative","cytokine","others"))
ggscatter(exprCPKgroup,x='corr',y='log10padj',color='CtyLabels',palette=c(get_palette(palette='npg',3),"grey"),label='expr',repel=T,label.select=cytokinePick$expr,
  font.label = c(12, "bold"), size = 1.5,
  xlab="correlation between gene expression and (IgG1+IgG3)%" ) + geom_vline( xintercept = 0.1, linetype = "dashed") +
  geom_vline( xintercept = -0.1, linetype = "dashed") + my_theme
dev.off()

##===========check some genes========
smad = subset(exprCPKgroup, grepl("SMAD", expr))
runx = subset(exprCPKgroup, grepl("RUNX", expr))
tgfb = subset(exprCPKgroup, expr == "TGFB1")
pickGene = subset( exprCPK,  grepl("SMAD", expr) | grepl("RUNX", expr) | expr == "TGFB1")
pdf("Fig_TGFB_SMAD_RUNX_IgG13_fraction_in_each_cancer_type.pdf", width = 15, height = 8)
print( ggscatter( pickGene, x = "cor", y = "pvalue", facet.by = "expr", scales = "free", legend = "none", color = "Disease", label = "Disease", font.label = 10, repel = T, title = i) )
dev.off()

##========check cytokines in each cancer type============
cytokineGene = subset(exprCPK, expr %in% cytokine$V1)
pdf("Fig_heatmap_cytokine_TPM_IgG13_fraction.pdf",width=18,height=25)
corPlot(cytokineGene,'Disease','expr',rho='cor',"padj",max_value = max(cytokineGene$cor,na.rm=T),min_value = min(cytokineGene$cor,na.rm=T), box_size = 3, bar_anno = c(min(cytokineGene$cor,na.rm=T), 0 , max(cytokineGene$cor,na.rm=T)))
## remove some mucosal cancer types
mucosalCancer = c("LUSC","LUAD","CESC","COAD","READ","ESCA","STAD","UCEC")
focusCancer = c("BLCA","OV","SKCM")
subdata1 = subset( cytokineGene, ! Disease %in% mucosalCancer )
corPlot(subdata1,'Disease','expr',rho='cor',"padj",max_value = max(subdata1$cor,na.rm=T),min_value = min(subdata1$cor,na.rm=T), box_size = 3, bar_anno = c(min(subdata1$cor,na.rm=T), 0 , max(subdata1$cor,na.rm=T)))
## focus on BLCA, OV and SKCM
subdata2 = subset( cytokineGene,  Disease %in% focusCancer )
corPlot(subdata2,'Disease','expr',rho='cor',"padj",max_value = max(subdata2$cor,na.rm=T),min_value = min(subdata2$cor,na.rm=T), box_size = 3, bar_anno = c(min(subdata2$cor,na.rm=T), 0 , max(subdata2$cor,na.rm=T)))
dev.off()

##======top7 gene in each cancer type=======
top7 = exprCPK %>% group_by(Disease) %>% top_n(7,cor)
top7expr = subset(exprCPK, expr %in% top7$expr)
pdf("Fig_heatmap_top7eachCancer_TPM_IgG13_fraction.pdf",width=18,height=25)
corPlot(top7expr,'Disease','expr',rho='cor',"padj",max_value = max(cytokineGene$cor,na.rm=T),min_value = min(cytokineGene$cor,na.rm=T), box_size = 3, bar_anno = c(min(cytokineGene$cor,na.rm=T), 0 , max(cytokineGene$cor,na.rm=T)))
# remove mucosal cancer types
subdata1 = subset( top7expr, ! Disease %in% mucosalCancer )
corPlot(subdata1,'Disease','expr',rho='cor',"padj",max_value = max(subdata1$cor,na.rm=T),min_value = min(subdata1$cor,na.rm=T), box_size = 3, bar_anno = c(min(subdata1$cor,na.rm=T), 0 , max(subdata1$cor,na.rm=T)))
## focus on BLCA, OV and SKCM
subdata2 = subset( top7expr,  Disease %in% focusCancer )
corPlot(subdata2,'Disease','expr',rho='cor',"padj",max_value = max(subdata2$cor,na.rm=T),min_value = min(subdata2$cor,na.rm=T), box_size = 3, bar_anno = c(min(subdata2$cor,na.rm=T), 0 , max(subdata2$cor,na.rm=T)))
dev.off()

##======top200 gene in each cancer type=======
exprCPKsig = subset( exprCPK, pvalue < 0.05)
top200 = exprCPKsig %>% group_by(Disease) %>% arrange(-cor) %>% top_n(200,cor)
write.table( top200, "top200_IgG13_fraction_genes_for_GO.txt", sep="\t", quote=F, row.names=F, col.names=T)
