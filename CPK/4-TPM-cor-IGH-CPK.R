library(dplyr)
library(data.table)
library(tidyr)
library(survcomp)
library(ggpubr)
library(ggplot2)
setwd('/liulab/zouyang/TCGA/IGH/6_trust4_clustering/TPM_corr_CPK')
source("/liulab/zouyang/TCGA/IGH/NewCgene/SHMratio_res/heatmap.r")
my_theme = theme(legend.title = element_blank(),
                 legend.text = element_text(size=17, face = "bold"),
                 legend.key.width=unit(1, "lines"),
                 plot.title = element_text(size=10, face="bold", hjust=0.5),
                 #plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),
                 axis.text.x = element_text(size=17, face="bold"),
                 axis.title.x=element_text(size=17, face="bold"),
                 axis.text.y = element_text(size = 17, face="bold"),axis.title.y = element_text(size = 17,face="bold"),
                 panel.grid.major = element_blank())

mucosalCancer = c("LUSC","LUAD","CESC","COAD","READ","ESCA","STAD","UCEC")
## correlation between CPK and gene expression
corr_path = './'
exprCPK = NULL
for(f in list.files(corr_path, pattern='corr_TPM_and_CPKcomplete_spearman.*.Rdata')) {
    load(paste0(corr_path, f))
    if( length(exprCPK.pcor)>1 ){
      cancerType = gsub('corr_TPM_and_CPKcomplete_spearman.','',f)
      cancerType = gsub('.old.Rdata','',cancerType)
      if( cancerType %in% mucosalCancer ){ next }
      exprCPK = rbind(exprCPK, data.frame(exprCPK.pcor,Disease=cancerType,expr=rownames(exprCPK.pcor)))
    }
}
save(exprCPK,file="./AllCorr_TPM_and_CPKcomplete_spearman.IGH.old.Rdata")
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
save(exprCPKpvalSort,file='./geneOrder_CPKcomplete.IGH.old.Rdata')
write.table(na.omit(exprCPKpvalSort[,c('expr','medianCor')]),"gene_corr_CPK_IGH.rnk",sep='\t',row.names=F,quote=F)

## rank gene by CPK
top10CPK = exprCPKpvalSort$expr[1:cutoff]
top10 = as.character(exprCPKpvalSort$expr[1:20])
write.table(as.character(exprCPKpvalSort$expr[1:200]),"top200_genes_in_CPKcomplete.IGH.old.txt",sep=",",quote=F,row.names=F,col.names=F)
exprCPKgroup <- exprCPK %>% group_by(expr) %>% summarise(corr=median(cor,na.rm=T),pvalue=combine.test(pvalue,na.rm=T))
exprCPKgroup$padj = p.adjust(exprCPKgroup$pvalue, method='BH')
exprCPKgroup = exprCPKgroup[order(exprCPKgroup$padj),]
exprCPKgroup$log10padj = -log10(exprCPKgroup$padj)
exprCPKgroup$labels = sapply(exprCPKgroup$expr, function(x) {
  if(x %in% top10){
    return('top')
  }else{
    return('others')
  }
})
pdf('Fig_rank_gene_by_correlation_TPM_CPKcomplete.IGH.pdf')
exprCPKgroup$labels <-factor(exprCPKgroup$labels, levels=c("top","others"))
ggscatter(exprCPKgroup,x='corr',y='log10padj',color='labels',label='expr',repel=T,label.select=c(top10),
          font.label = c(12, "bold"), size = 1.5,
          xlab="correlation between gene expression and IGH diversity",ylab="-log10(Padj)",palette=c(get_palette(palette='npg',1),"grey"))+ my_theme
dev.off()

## heatmap of top 10 genes
exprCPK.sub = subset(exprCPK, expr %in% top10CPK)
exprCPK.sub$expr = factor(exprCPK.sub$expr, levels=top10CPK)
maxVal = round(max(exprCPK$cor,na.rm=T),2)
minVal = round(min(exprCPK$cor,na.rm=T),2)
exprCPK$cor[ exprCPK$cor > maxVal] = maxVal
exprCPK$cor[ exprCPK$cor < minVal] = minVal

p1<-corPlot(exprCPK.sub,'Disease','expr',rho='cor',"pvalue",max_value = maxVal,min_value = minVal, box_size = 5.5, bar_anno = c(minVal, 0 , maxVal))
ggsave(p1,file="./Fig_heatmap_TPM_CPKcomplete.IGH.pdf",width=7,height=4,dpi=1000)
