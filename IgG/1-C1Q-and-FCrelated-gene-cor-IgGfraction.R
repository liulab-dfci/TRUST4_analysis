library(ppcor)
library(reshape2)
library(ggpubr)
#setwd("/liulab/zouyang/TCGA/IGH/6_trust4_clustering/TPM_corr_CPyK")
setwd("/liulab/lsong/projects/TRUST4_Ouyang_analysis/Fig5")
load("/liulab/zouyang/TCGA/IGH/1_SHM/sampleInfo_and_cancer_type.Rdata")
load("/liulab/zouyang/TCGA/IGH/1_SHM/Ig_mat.Rdata")
Ig.matProp = prop.table(as.matrix(Ig.mat[,3:11]),margin=1)
Ig.matProp = as.data.frame(Ig.matProp)
rownames(Ig.matProp) = Ig.mat$TCGA_id
Ig.mat = Ig.matProp
Ig.mat$IgG1_plus_IgG3 = Ig.mat[,'IGHG1']+Ig.mat[,'IGHG3']
Ig.mat$IGHA = Ig.mat[,'IGHA1']+Ig.mat[,'IGHA2']
Ig.mat$IGHG = Ig.mat[,'IGHG1'] + Ig.mat[,'IGHG2'] + Ig.mat[,'IGHG3'] + Ig.mat[,'IGHG4']

#============ annotation information==============
## uniform digit of sample id
getID<-function(sId,digits=3){
  tmp=sapply(sId,function(x){
    return(
      paste0(
        strsplit(x,split = '\\-')[[1]][1:digits],collapse='-')
    )})
  return(tmp)
}

# Trust4's cancer type
tumorPurity=get(load('/liulab/lsong/data/bcr/BCR_NG_TRUST/data/TRUST3/tumorPurity.Rdata'))
exprDat = get(load('/liulab/zouyang/TCGA/QuantFile/TPM_of_TCGA.proteinCoding.Rdata'))

#================calculate correlation======================
#testGene = c("FCGR1A","FCGR2A","FCGR2B","FCGR2C","FCGR3A","FCGR3B",'MS4A2', "FCER2", 'FCAR', "FCGRT", "FCAMR", "PIGR")
testGene = c("C2","C3","C4A","C4B","C5","C6","C7",'C8A', "C8B", 'C8G', "C9")
tmp.expr = exprDat[testGene,]
MHCII = apply(exprDat[grep('^HLA-D',rownames(exprDat)),],2,mean)
MHCI = apply(exprDat[grep('^HLA-[ABC]',rownames(exprDat)),],2,mean)
C1Q = apply(exprDat[c("C1QA", "C1QB", "C1QC"),],2,mean)
tmp.expr = rbind( tmp.expr, MHCII = MHCII, MHCI = MHCI, C1Q = C1Q)
#testGene = c("C1Q","FCGR1A","FCGR2A","FCGR2B","FCGR2C","FCGR3A","FCGR3B",'MS4A2', "FCER2", 'FCAR', "FCGRT", "FCAMR", "PIGR")

tmp.cor.List=c()
cancerOrder = sort(unique(cancer.types))
for(cc in cancerOrder){
  #cc = "BRCA"
  print(cc)
	tmp.ss = Reduce( intersect, list( names(cancer.types[which(cancer.types==cc)]), colnames(tmp.expr), rownames(Ig.mat)))
	tmp.cor.mat = matrix(NA,nrow= ncol(Ig.mat),ncol=nrow(tmp.expr))
	rownames(tmp.cor.mat) = colnames(Ig.mat)
	colnames(tmp.cor.mat) = rownames(tmp.expr)
	tmp.p.mat=tmp.cor.mat
	for(ii in 1:ncol(Ig.mat)){
    #ii = 12
		for(jj in 1:nrow(tmp.expr)){
      #jj = 10 # FCGRT
			tmp.dd0= as.data.frame( cbind( t(tmp.expr[jj, tmp.ss]),Ig.mat[tmp.ss, ii]) )
      tmp.dd0$tumorPurity = tumor.purity[getID(tmp.ss,4)]
			tmp.dd0=na.omit(tmp.dd0)
			if(nrow(tmp.dd0) < 30)next
			tmp=pcor.test(tmp.dd0[,1],tmp.dd0[,2],tmp.dd0[,3],method='s')
			tmp.cor.mat[ii,jj]=tmp$estimate
			tmp.p.mat[ii,jj]=tmp$p.value
		}
	}
	tmp.cor.List=c(tmp.cor.List,list(list(COR=tmp.cor.mat,Pval=tmp.p.mat)))
}
names(tmp.cor.List)=cancerOrder
save(testGene, tmp.cor.List,file = "CorList_of_C1Q_Genes_and_IgFraction.Rdata")

#===============plot heatmap====================

source('/homes/zouyang/Immunotherapy/src/heatmap_size.R')
testIsotype=c("IGHE", "IGHD", 'IGHA','IGHA2','IGHA1','IGHG4','IGHG3','IGHG2','IGHG1','(IgG1+IgG3)%','IGHG',"IGHM")
testIsotype=c( 'IGHA2', "IGHE", 'IGHG4','IGHG2','IGHA1','IGHG1','IGHG3',"IGHD","IGHM")
testIsotype=c('IGHA2','IGHG4','IGHG2','IGHA1','IGHG1','IGHG3',"IGHM")
# all receptor genes
#testGene = c("FCGR1A","FCGR2A","FCGR2B","FCGR3A","FCGR3B",'MS4A2', "FCER2", 'FCAR', "FCGRT", "FCAMR", "PIGR", "C1Q")#no expression for "FCGR2C"
testGene = c("C2","C3","C4A","C4B","C5","C6","C7",'C8A', "C8B", 'C8G', "C9")
pdf("Fig_C2-C9gene_correlated_with_IgFraction.pdf", width=12, height=15)
for(gene in testGene){
    #print( gene )
    corGene.mat=c()
    corGene.pmat=c()
    for(cc in names(tmp.cor.List)){
        #print(cc)
        corGene.mat=cbind(corGene.mat,tmp.cor.List[[cc]][[1]][,gene])
        corGene.pmat=cbind(corGene.pmat,tmp.cor.List[[cc]][[2]][,gene])
    }
    colnames(corGene.mat)=colnames(corGene.pmat)=names(tmp.cor.List)

    corGene.mat = melt(corGene.mat)
    corGene.pmat = melt(corGene.pmat)
    corGene.mat$pvalue = corGene.pmat$value
    corGene.mat=na.omit(corGene.mat)
    if( nrow(corGene.mat)==0 ){next}

    ###Check for Isoytoe###
    test.mat = subset(corGene.mat, Var1  %in% testIsotype )
    test.mat$Padj = p.adjust(test.mat$pvalue, method = 'BH')
    test.mat$Var1 = sapply(as.character( test.mat$Var1 ), function(x){ if(x=="IgG1_plus_IgG3"){return( "(IgG1+IgG3)%" )}else{return(x)}})
    test.mat$Var1 = factor(test.mat$Var1,levels = testIsotype)

    ### gene name alias
    if( gene == "FCGR1A"){
      geneName = "CD64"
    }else if( gene == "FCGR2A" ){
      geneName = "CD32A"
    }else if( gene == "FCGR2B" ){
      geneName = "CD32B"
    }else if( gene == "FCGR2C" ){
      geneName = "CD32C"
    }else if( gene == "FCGR3A" ){
      geneName = "CD16a"
    }else if( gene == "FCGR3B" ){
      geneName = "CD16b"
    }else if( gene == "FCGRT" ){
      geneName = "FcRN"
    }else if( gene == "MS4A2" ){
      geneName = "FCER"
    }else if( gene == "CD274" ){
      geneName = "PDL1"
    }else if( gene == "PDCD1" ){
      geneName = "PD1"
    }else{
      geneName = gene
    }

    minVal = round(min(test.mat$value),2)
    maxVal = round(max(test.mat$value),2)
    print( geneName )

    test.mat$value[ test.mat$value>maxVal ] = maxVal
    test.mat$value[ test.mat$value<minVal ] = minVal
    print( corPlot(dat = test.mat,x_axis = 'Var2',y_axis = 'Var1',rho = 'value',pvalue = 'Padj',box_size = 5,
           max_value= maxVal, min_value= minVal,titleName = geneName,
           bar_anno  = c( minVal, 0, maxVal) ) )

}
dev.off()
quit()

#===============FcRN and survival in two class  ============
fcrnExpr = unlist(exprDat['FCGRT',])
result<-generateData(igRatio = fcrnExpr,
                     clinDataPath = '/liulab/lsong/data/bcr/BCR_NG_TRUST/data/TRUST3/clin.Rdata',
                     cancerTypePath='/liulab/zouyang/TCGA/IGH/1_SHM/sampleInfo_and_cancer_type.Rdata',
                     tumorPurityPath='/liulab/lsong/data/bcr/BCR_NG_TRUST/data/TRUST3/tumorPurity.Rdata')

pdf('figure_survival_FcRN_FPKM.pdf',width=15,height=20)
## reproduce figure6d
opar = par(par(no.readonly = T))
par(mfrow=c(7,5))
for(i in unique(result$cancer)){
  res_brca <- subset(result,cancer==i)
  figure6d(res_brca,title=i)
}
par(opar)
dev.off()

#===============FcRN and survival in four class  ============
source("func_survival_4class.R")
## IgG1
igg1 = Ig.mat[, "IGHG1"]
names(igg1) = rownames(Ig.mat)
result<-generate2Data(shm = fcrnExpr,
                     clinDataPath = '/liulab/lsong/data/bcr/BCR_NG_TRUST/data/TRUST3/clin.Rdata',
                     igFraction = igg1,
                     cancerTypePath='/liulab/zouyang/TCGA/IGH/1_SHM/sampleInfo_and_cancer_type.Rdata',
                     tumorPurityPath='/liulab/lsong/data/bcr/BCR_NG_TRUST/data/TRUST3/tumorPurity.Rdata')

pdf('figure_survival_FcRN_FPKM_and_IgG1Fraction.pdf',width=15,height=20)
## reproduce figure6d
opar = par(par(no.readonly = T))
par(mfrow=c(7,5))
for(i in unique(result$cancer)){
  res_brca <- subset(result,cancer==i)
  figure5c(res_brca,title=i)
}
par(opar)
dev.off()
