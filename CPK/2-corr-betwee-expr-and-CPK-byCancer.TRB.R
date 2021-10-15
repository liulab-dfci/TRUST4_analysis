p = c('ppcor','reshape2')
for(el in p){
  if (!is.element(el, installed.packages()[,1]))install.packages(el, dep=TRUE)
  suppressWarnings(suppressMessages(invisible(require(el, character.only=TRUE))))
}

# clean env
rm(list = ls())
########################Function#########################################
## uniform digit of sample id
getID<-function(sId,digits=3){
  tmp=sapply(sId,function(x){
    return(
      paste0(
        strsplit(as.character(x),split = '\\-')[[1]][1:digits],collapse='-')
    )})
  return(tmp)
}

getGeneExprs<-function(exprDat,genelist){
  sp.Gene = c('MHCII','MHCI','CD16a')
  com.Expr = exprDat[setdiff(genelist,sp.Gene),,drop=F]
  exist.sp.Gene = intersect(genelist,sp.Gene)
  #### Return corresponding gene expression if only extract based on gene name
  if(length(exist.sp.Gene)<1)return(com.Expr)

  #### Extract special gene expression
  sp.Expr=c()
  for (gene in exist.sp.Gene) {
      if(gene == 'MHCII'){
        tmp.Expr = apply(exprDat[grep('^HLA-D',rownames(exprDat)),],2,mean)
      } else if(gene == 'MHCI'){
        tmp.Expr = apply(exprDat[grep('^HLA-[ABC]',rownames(exprDat)),],2,mean)
      } else {
        tmp.Expr = exprDat['FCGR3A',]
      }
    sp.Expr=rbind(sp.Expr,tmp.Expr)
  }
  rb.sp_com.Expr = rbind(sp.Expr,com.Expr)
  rownames(rb.sp_com.Expr)[1:length(exist.sp.Gene)] = exist.sp.Gene
  return(rb.sp_com.Expr)
}


combineDat <-function(swithDat,cancerType,tumorPurity,exprDat,testGene){
  ### modify name
  names(swithDat) = getID(names(swithDat),4)
  names(cancerType)=getID(names(cancerType),4)
  exist.Cancer = cancerType[names(swithDat)]
  colnames(exprDat)= getID(colnames(exprDat),4)

  ## log2 transform expression valu
  exprDat= log2(exprDat+1)

  ## get intersect sample
  overlapSample = intersect(colnames(exprDat),names(exist.Cancer))
  swithDat=swithDat[overlapSample]
  tumorPurity = tumorPurity[overlapSample]
  exist.Cancer =exist.Cancer[overlapSample]
  ## get needed expression data
  testExprs = getGeneExprs(exprDat = exprDat,testGene)
  exprSwith = as.data.frame(t(rbind(testExprs[,overlapSample],
                                            Purity = tumorPurity,
                                            Ig1_Ig3=swithDat,
                                            cancer=as.character(exist.Cancer))),stringsAsFactors = F)
  return(exprSwith)

}

getExprSwithPcor <-function(exprSwith,genelist,selectCancerType){

  sub.exprSwtich  = exprSwith[exprSwith$cancer == selectCancerType,-match('cancer',colnames(exprSwith))]
  #if(dim(sub.exprSwtich)[1]<30 || length(unique(sub.exprSwtich$Ig1_Ig3))==1){return(NA)}
  sub.exprSwtich = apply(sub.exprSwtich, 2, as.numeric)
  corMatrix = matrix(NA, nrow=length(genelist), ncol = 2,
                   dimnames=list(genelist,c('cor','pvalue')))
  for(GENE in genelist){
      print(GENE)
      exprSwitchGene = na.omit(sub.exprSwtich[,c('Ig1_Ig3','Purity',GENE)])
      if( nrow(exprSwitchGene)>=30 && length(unique(exprSwitchGene[,'Ig1_Ig3']))>1 && length(unique(exprSwitchGene[,GENE]))>1 && var(exprSwitchGene[,GENE])>.Machine$double.eps){
        corValue=pcor(exprSwitchGene[,c('Ig1_Ig3','Purity',GENE)],method='spearman')
        corMatrix[GENE,'cor']=corValue$estimate['Ig1_Ig3',GENE] # cor.cor
        corMatrix[GENE,'pvalue']=corValue$p.value['Ig1_Ig3',GENE] #pvalue.pvalue
      }
  }
  return(corMatrix)
}
########################Function#########################################


######### MAIN ##############
## INPUT
setwd('/liulab/zouyang/TCGA/IGH/6_trust4_clustering/TPM_corr_CPK')
print('start...')
load("/liulab/zouyang/TCGA/IGH/1_SHM/metrics_of_TRB_and_IGH.rmSingleton.Rdata")
swithDat <- igh_trb$TRB.CPK
names(swithDat) <- igh_trb$ID

# Trust4's cancer type
load('/liulab/zouyang/TCGA/IGH/1_SHM/sampleInfo_and_cancer_type.Rdata')
tumorPurity=get(load('/liulab/lsong/data/bcr/BCR_NG_TRUST/data/TRUST3/tumorPurity.Rdata'))
#exprDat = readRDS('/liulab/lsong/data/bcr/BCR_NG_TRUST/data/TRUST3/PanCancerTMP.rds')
exprDat = get(load('/liulab/zouyang/TCGA/QuantFile/TPM_of_TCGA.proteinCoding.Rdata'))
######### Gene you want to test #########
testGene <- rownames(exprDat)
exprSwith = combineDat(swithDat,cancerType,tumorPurity,exprDat,testGene)
args = commandArgs(trailingOnly = TRUE)
selectCT <- args[1]
print('running...')
exprCPK.pcor = getExprSwithPcor(exprSwith = exprSwith,genelist =testGene,selectCancerType = selectCT)
save(exprCPK.pcor,file=paste0("corr_TPM_and_TRB_CPK_rmSingleton_spearman.",selectCT,".Rdata"))
