setwd("~/Dropbox (Partners HealthCare)/Immunotherapy/inter_cohorts/")
load("CPK_allData.Rdata")
colnames(allIgProp)[1] = "Patient"
allData$normUniqCDR3 = allData$uniqCDR3/allData$totalReads
allData$normT4CDR3 = allData$libsize/allData$totalReads
allIgProp$IgG1_plus_IgG3 = allIgProp$IGHG1 + allIgProp$IGHG3
load("../htseq_vdj_from_David/TPM_of_CTLA4_and_PDL1.Rdata")
load("../htseq_vdj_from_David/immuneDeconv/All_FC_of_immuneDeconv_in_response.Rdata")
metaData = Reduce( function(x,y){ merge(x,y,by=c("Patient","cohort"))}, list( allData, allIgProp, allTPMdata, allImmuneDeconv))
unique(metaData$cohort)
colnames(metaData) = sapply(colnames(metaData), function(x){
  if( x == "CPK"){
    return("Diversity")
  }else if(x=="normUniqCDR3"){
    return("Richness")
  }else if(x=="shm"){
    return("SHM")
  }else if(x=="normT4CDR3"){
    return("TRUST4_IGHabundance")
  }else if(x=="normUniqCluster"){
    return("uniqCluster.norm")
  }else if(x == "CD274"){
    return("PDL1")
  }else if(x == "PDCD1"){
    return("PD1")
  }else{
    return(x)
  }
})
testMetric = c("IgG1_plus_IgG3","Richness","Diversity","TRUST4_IGHabundance","B.cell_MCP_COUNTER","T.cell.CD8._MCP_COUNTER","PDL1","PD1","CTLA4")

#===========mixed effect linear regression model ==========
subMetaData = subset(metaData, !grepl("SQ|LN|bladder|kidney|lymph|ureter|Gide_PD1|Gide_ipiPD1",cohort) )
unique(subMetaData$cohort)
subMetaData$cancer = sapply( subMetaData$cohort, function(x){
      if( grepl("Mel", x) ){
        return("SKCM")
      }else if( grepl("mRCC",x)){
        return("KIRC")
      }else if( grepl("mGC",x)){
        return("STAD")
      }else if( grepl("mUC",x)){
        return("BLCA")
      }
})
subMetaData = subset( subMetaData, libsize > quantile(subMetaData$libsize,0.05))
preMetadata = subset( subMetaData, !grepl("On|EDT", cohort))
unique(preMetadata$cohort)
onMetadata = subset( subMetaData, grepl("On|EDT", cohort))
unique(onMetadata$cohort)

#==========whether covariates are correlated===========
library(corrplot)
pdf("Fig_cor_between_metrics.pdf")
col3 <- colorRampPalette(c("navy", "white", "firebrick3"))
preCor = cor( preMetadata[, testMetric], method = "spearman" , use = "complete.obs")
corrplot(preCor, method = "number", col = col3(50), tl.col = "black", title = "Pre")
onCor = cor( onMetadata[, testMetric], method = "spearman" , use = "complete.obs")
corrplot(onCor, method = "number", col = col3(50), tl.col = "black", title = "On")
dev.off()
#==============mixed effects logistic regression=========
library(lmerTest)
#install.packages("lme4")
#install.packages("glmerTest")
#=============pre data===============
quantileRange <- function(data,metric,cutoff = 0.05){
  #data = subMetaData
  #metric = "Diversity"
  lowB = quantile(data[,metric], cutoff)
  highB = quantile(data[,metric], 1-cutoff)
  data = data[which(data[,metric]>lowB & data[,metric] <highB), ]
  return(data)
}
mixedLogisticModel <- function( preMetadata){
  glmerPval = NULL
  full.glmer = glmer(Response.x ~ IgG1_plus_IgG3 +  (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IgG1_plus_IgG3"))

  full.glmer = glmer(Response.x ~ IGHG1 +  (1|cohort), data = preMetadata, family = binomial,control = glmerControl(optimizer = "bobyqa"))
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHG1"))

  full.glmer = glmer(Response.x ~ IGHG2 +  (1|cohort), data = preMetadata, family = binomial,control = glmerControl(optimizer = "bobyqa"))
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHG2"))

  full.glmer = glmer(Response.x ~ IGHG3 +  (1|cohort), data = preMetadata, family = binomial,control = glmerControl(optimizer = "bobyqa"))
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHG3"))

  full.glmer = glmer(Response.x ~ IGHG4 +  (1|cohort), data = preMetadata, family = binomial,control = glmerControl(optimizer = "bobyqa"))
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHG4"))

  full.glmer = glmer(Response.x ~ IGHA1 +  (1|cohort), data = preMetadata, family = binomial,control = glmerControl(optimizer = "bobyqa"))
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHA1"))

  full.glmer = glmer(Response.x ~ IGHA2 +  (1|cohort), data = preMetadata, family = binomial,control = glmerControl(optimizer = "bobyqa"))
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHA2"))

  full.glmer = glmer(Response.x ~ IGHD +  (1|cohort), data = preMetadata, family = binomial,control = glmerControl(optimizer = "bobyqa"))
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHD"))

  full.glmer = glmer(Response.x ~ IGHM +  (1|cohort), data = preMetadata, family = binomial,control = glmerControl(optimizer = "bobyqa"))
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHM"))

  full.glmer = glmer(Response.x ~ IGHE +  (1|cohort), data = preMetadata, family = binomial,control = glmerControl(optimizer = "bobyqa"))
  res = summary( full.glmer )
  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHE"))

  full.glmer = glmer(Response.x ~ log2(Richness) + (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "Richness"))

  full.glmer = glmer(Response.x ~ log2(Diversity) + (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "Diversity"))


  full.glmer = glmer(Response.x ~ log2(TRUST4_IGHabundance) +  (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "TRUST4_IGHabundance"))


  full.glmer = glmer(Response.x ~ B.cell_MCP_COUNTER +  (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "B.cell_MCP_COUNTER"))


  full.glmer = glmer(Response.x ~ T.cell.CD8._MCP_COUNTER +  (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "T.cell.CD8._MCP_COUNTER"))


  full.glmer = glmer(Response.x ~ B.cell_XCELL +  (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "B.cell_XCELL"))


  full.glmer = glmer(Response.x ~ log2(PDL1+1) +  (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "PDL1"))

  full.glmer = glmer(Response.x ~ log2(PD1+1) +  (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "PD1"))

  full.glmer = glmer(Response.x ~ log2(CTLA4+1) +  (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "CTLA4"))

  full.glmer = glmer(Response.x ~ log2(IGHC+1) +  (1|cohort), data = preMetadata, family = binomial)
  res = summary( full.glmer )

  glmerPval = rbind( glmerPval, data.frame( coef = res$coefficients[2,1], pval = res$coefficients[2,4], metric = "IGHC"))

  glmerPval = glmerPval[,c(3,1,2)]
  return(glmerPval)
}

# filter sample by IGH lib size
prePval = mixedLogisticModel( preMetadata = preMetadata)
onPval = mixedLogisticModel( preMetadata = onMetadata )
allPval = rbind( data.frame( prePval, condition = "Pre"),
                 data.frame( onPval, condition = "On"))
allPval = subset( allPval, metric %in% testMetric)
allPval$metric = sapply(as.character(allPval$metric), function(x){
  if(x=="IgG1_plus_IgG3"){
    return("(IgG1+IgG3)%")
  }else{
    return(x)
  }
})

allPval$padj = p.adjust( allPval$pval, method = "BH" )
maxVal= 1
minVal = -1
allPval$coefSign[ allPval$coef > 0 ] = 1
allPval$coefSign[ allPval$coef < 0 ] = -1
allPval$metric = factor( allPval$metric, levels = c("(IgG1+IgG3)%","Richness","Diversity","TRUST4_IGHabundance","B.cell_MCP_COUNTER","T.cell.CD8._MCP_COUNTER","PDL1","PD1","CTLA4"))
pdf("Fig_mixed_model_for_combinedCohorts.pdf")
corPlot( allPval, x_axis = "condition", y_axis = "metric", rho = "coefSign", pvalue = "padj", title = "",
         fdr_cutoff = 0.1, min_value = minVal, max_value = maxVal, mean_value = 0, bar_anno=c(minVal,0,maxVal), theme = mytheme_grid )
dev.off()
