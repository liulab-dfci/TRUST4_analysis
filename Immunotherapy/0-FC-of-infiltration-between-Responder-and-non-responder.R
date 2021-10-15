FCinfiltration <- function(infilPath, clinicalPath, cohortName){
  # Fold Change of infiltration between Responder and Non-responder
  #infilPath = "Mel_PD1_Riaz_Pre_ipiProg.immuneDeconv"
  #clinicalPath = "../clinical/Mel_PD1_Riaz_Pre_ipiProg_clinical.Rdata"
  #cohortName = "Mel_PD1_Riaz_Pre_ipiProg"
  infilData = read.table(infilPath, sep="\t",row.names=1,header=T, stringsAsFactor=F)
  clinicalData = get(load(clinicalPath))
  clinicalData = clinicalData[ !is.na(clinicalData$Response), ]
  # some cohorts's sample names are not consistent
  overlapSample = intersect(colnames(infilData), clinicalData$Patient)
  overlapId = "Patient"
  if( length(overlapSample)==0 ){
    overlapId = "Run"
  }
  clinicalData = clinicalData[ !is.na(clinicalData[, overlapId]), ]
  colAnno = data.frame(response = factor(clinicalData$Response, levels = c(1,0)))
  rownames(colAnno) = clinicalData[, overlapId]

  R = subset(clinicalData, Response==1, select=overlapId)
  NR = subset(clinicalData, Response==0, select=overlapId)
  R.infil = infilData[,colnames(infilData) %in% R[,overlapId]]
  NR.infil = infilData[,colnames(infilData) %in% NR[,overlapId]]

  # heatmap of MCP_COUNTER
  subdata = cbind(R.infil, NR.infil)
  subdata = subset( subdata, grepl("MCP_COUNTER",rownames(infilData)))
  #pheatmap(na.omit( subdata ), scale="row", annotation_col =colAnno, cluster_cols=F, gaps_col = ncol(R.infil))
  #pheatmap(subdata, scale="row", annotation_col =colAnno, gaps_col = ncol(R.infil))

  # fold change in response
  FCmat = Reduce(rbind, lapply(rownames(infilData), function(x){
        Rvalue =  as.numeric(R.infil[x,])
        NRvalue = as.numeric(NR.infil[x,])
        tmp.FC = median( Rvalue ,na.rm=T)/median( NRvalue, na.rm=T)
        #tmp.P = NA
        #if( length(Rvalue)>2 & length(NRvalue)>2){
          tmp.P = wilcox.test(Rvalue,NRvalue)$p.value
      #  }
        return(data.frame(cellType = x, cohort = cohortName, FC = tmp.FC, Pvalue = tmp.P))
  }) )
  return( FCmat )
}

combneInfiltration <- function(infilPath, clinicalPath, cohortName){
  # Fold Change of infiltration between Responder and Non-responder
  #infilPath = "Mel_PD1_Riaz_Pre_ipiProg.immuneDeconv"
  #clinicalPath = "../clinical/Mel_PD1_Riaz_Pre_ipiProg_clinical.Rdata"
  #cohortName = "Mel_PD1_Riaz_Pre_ipiProg"
  infilData = read.table(infilPath, sep="\t",row.names=1,header=T, stringsAsFactor=F)
  clinicalData = get(load(clinicalPath))
  clinicalData = clinicalData[ !is.na(clinicalData$Response), ]
  # some cohorts's sample names are not consistent
  overlapSample = intersect(colnames(infilData), clinicalData$Patient)
  overlapId = "Patient"
  if( length(overlapSample)==0 ){
    overlapId = "Run"
    overlapSample = intersect(colnames(infilData), clinicalData$Run)
  }
  clinicalData = clinicalData[ !is.na(clinicalData[, overlapId]), ]
  rownames(clinicalData) = clinicalData[, overlapId]
  clinicalData = clinicalData[ overlapSample, ]
  colAnno = data.frame(response = factor(clinicalData$Response, levels = c(1,0)))
  rownames(colAnno) = clinicalData[, overlapId]
  pid = clinicalData
  combineData = data.frame( t(infilData[ grepl("MCP_COUNTER|XCELL",rownames(infilData) ), overlapSample]), Response = colAnno[ overlapSample, ], Patient = clinicalData$Patient, cohort = cohortName )

  return( combineData )
}

#------------analysis----------------------
#-----correlation between immune cell type and IGH CPK-------
setwd("/homes/zouyang/Immunotherapy/infiltration")
source("/homes/zouyang/Immunotherapy/src/heatmap_size.R")
library(pheatmap)

clinicalDir = "/homes/zouyang/Immunotherapy/clinical/"
infilDir = "/homes/zouyang/Immunotherapy/infiltration/"

cohorts = list.files(path = infilDir ,pattern = ".immuneDeconv$")
cohorts
allImmuneDeconv = Reduce(rbind, lapply(cohorts, function(x){
        #x = cohorts[6]
        cohort = gsub(pattern = ".immuneDeconv",replacement = "",x)
        print(cohort)
        clinicalFile = paste0(clinicalDir, cohort, "_clinical.Rdata")
        combData = combneInfiltration( infilPath = x, clinicalPath = clinicalFile, cohortName = cohort)
        print( dim(combData) )
        return( combData )
}))

allFCdata = Reduce(rbind, lapply(cohorts, function(x){
        #x = cohorts[6]
        cohort = gsub(pattern = ".immuneDeconv",replacement = "",x)
        print(cohort)
        clinicalFile = paste0(clinicalDir, cohort, "_clinical.Rdata")
        FCdata = FCinfiltration( infilPath = x, clinicalPath = clinicalFile, cohortName = cohort)
}))
allFCdata$log2FC = log2(allFCdata$FC)
allFCdata$Padj = p.adjust(allFCdata$Pvalue, method="BH")
#effectCohorts = c("Mel_ICB_Gide_ipiPD1_PRE_LN","Mel_ICB_Gide_PD1_PRE_SQ","Mel_ICB_Gide_PD1_PRE_LN","Mel_PD1_Riaz_On_ipiProg","Mel_PD1_Riaz_Pre_ipiProg","mRCC_PDL1_McDermott_Kidney")
effectCohorts = c("Mel_ICB_Gide_ipiPD1_PRE_LN","Mel_ICB_Gide_PD1_PRE_SQ","Mel_ICB_Gide_PD1_EDT_SQ","Mel_ICB_Gide_PD1_PRE_LN","Mel_PD1_Riaz_On_ipiProg","Mel_PD1_Riaz_Pre_ipiProg","Mel_PD1_Riaz_On_ipiNaive","Mel_PD1_Riaz_Pre_ipiNaive","mRCC_PDL1_McDermott_Kidney")
nonEffectCohorts = setdiff( as.character( unique(allFCdata$cohort) ), effectCohorts)
cohortOrder = c(effectCohorts, nonEffectCohorts)
cohortOrder
allFCdata$cohort = factor(allFCdata$cohort, levels=cohortOrder)
save(allFCdata, allImmuneDeconv, file = "All_FC_of_immuneDeconv_in_response.Rdata")
