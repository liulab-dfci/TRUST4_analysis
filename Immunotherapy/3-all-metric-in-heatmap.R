setwd('~/Dropbox (Partners HealthCare)/Immunotherapy/inter_cohorts/')
source("~/Dropbox (Partners HealthCare)/Immunotherapy/src/heatmap_size.R")
library(dplyr)
library(tidyr)
#==========rename cohorts===========
renameCohort <- function(x, cohortNum = NULL){
  x = gsub("^Mel","SKCM",x)
  x = gsub("^mRCC","KIRC",x)
  x = gsub("^mUC","BLCA",x)
  x = gsub("^mGC","STAD",x)
  if( !is.null(cohortNum) ){
    x = cohortName[ x ]
  }
  return(x)
}
scale_rows <- function(x){
  m = mean(x, na.rm = T)
  s = sd(x, na.rm = T)
  return((x - m) / s)
}

cauchyTest <- function(x){
  weight = 1/length(x)
  t = sum( weight*tan((0.5-x)*pi) )
  p = 1/2 - (atan(t))/pi
  return( p )
}

getFC <- function(metadata, metricName){
  cohortName = unique(metadata$cohort)
  tmp.FC = Reduce( rbind, lapply( cohortName, function(x){
    Res = subset(metadata, cohort==x & responseLabel == "responder")
    nonRes = subset(metadata, cohort==x & responseLabel == "non-responder")
    FC = median(Res[,metricName], na.rm=T)/median(nonRes[,metricName], na.rm=T)
    Pvalue = NA
    if( length(na.omit(Res[,metricName])) > 0 & length(na.omit(nonRes[,metricName]))>0){
      Pvalue = wilcox.test(na.omit(Res[,metricName]), na.omit(nonRes[,metricName]))$p.value
    }
    return( data.frame(cohort = x, metric = metricName, FC = FC, Pvalue = Pvalue) )
  }))

}
#===============IGH analysis=============
load("CPK_allData.Rdata")
load("SHM_of_cohorts.IGH.filterLowFreqCDR3.Rdata")
load("CSR_of_cohorts.IGH.Rdata")
load("cohort_color_for_plot.Rdata")
load("CohortT4Order_and_name_for_plot.Rdata")
load("../htseq_vdj_from_David/TPM_of_CTLA4_and_PDL1.Rdata")
##======sample number in each cohorts==========
allData$cohort = sapply(allData$cohort, renameCohort)
sampleNum = as.data.frame( table(allData$cohort, allData$responseLabel) )
pdf("Fig_patients_in_each_cohorts.pdf")
#sampleNum$Var1 = factor(sampleNum$Var1, levels = cohortT4Order)
ggbarplot(sampleNum, x = "Var1", y = "Freq", fill = "Var2", color = NA, label = TRUE, ylab = "number of patients")+my_theme_nocolor
sampleNum.sub = subset(sampleNum, !grepl("SQ|LN|bladder|kidney|lymph|ureter",Var1) )
sampleNum.sub = subset(sampleNum.sub, !grepl("Gide_PD1|Gide_ipiPD1",Var1) )
sampleNum.sub = subset(sampleNum.sub, !grepl("Riaz_Pre$|Riaz_On$",Var1) )
ggbarplot(sampleNum.sub, x = "Var1", y = "Freq", fill = "Var2", color = NA, label = TRUE, ylab = "number of patients")+my_theme_nocolor
sampleTotalNum = table(allData$cohort,allData$responseLabel)

cohortName = paste0(rownames(sampleTotalNum),"(",sampleTotalNum[,"responder"],",",sampleTotalNum[,"non-responder"],")")
names(cohortName) = rownames(sampleTotalNum)
cohortName = gsub("_Pre|_On|_PRE|_EDT|_all|Kidney","",cohortName)
cohortName = gsub("Mariathasan", "Mari", cohortName)
cohortName = gsub("McDermott", "McDe", cohortName)
cohortName = gsub("VanAllen", "VanA", cohortName)
#  cohortName to cohort
name2cohort = names(cohortName)
names(name2cohort) = cohortName

sampleNum.sub$cohort = sapply(sampleNum.sub$Var1, function(x){ renameCohort(x, cohortNum = cohortName) })
sampleNum.sub$cohort = factor(sampleNum.sub$cohort, levels = cohortT4Order)
colnames(sampleNum.sub)[2] = "response"
ggbarplot(sampleNum.sub, x = "cohort", y = "Freq", fill = "response", color = NA, label = TRUE, ylab = "number of patients")+my_theme_nocolor
dev.off()
save( renameCohort, name2cohort, cohortName, cohortT4Order, file = "function_rename_cohort.Rdata")

allData = subset(allData, uniqCDR3>0)
allData$normUniqCDR3 = allData$uniqCDR3/allData$totalReads
allData$normT4CDR3 = allData$libsize/allData$totalReads
cpkFC = Reduce( rbind, lapply(unique(allData$cohort), function(x){
          #x = "Mel_CTLA4_VanAllen"
          Res = subset(allData, cohort == x & responseLabel == "responder")
          nonRes = subset(allData, cohort == x & responseLabel == "non-responder")
          metric = c("entropy","evenness","CPK","normUniqCDR3","normT4CDR3","t2freq")
          tmp.FC = Reduce( rbind, lapply(metric, function(y){
                #y = "entropy"
                FC = median(Res[,y], na.rm=T)/median(nonRes[,y], na.rm=T)
                Pvalue = NA
                if( length(na.omit(Res[,y])) > 0 & length(na.omit(nonRes[,y]))>0){
                  Pvalue = wilcox.test(na.omit(Res[,y]), na.omit(nonRes[,y]))$p.value
                }
                return( data.frame(metric = y, FC = FC, Pvalue = Pvalue) )
          }))
          deltaEvenness = median(Res[, "evenness"], na.rm=T) - median(nonRes[, "evenness"], na.rm=T)
          Pvalue = wilcox.test(na.omit(Res[,"evenness"]), na.omit(nonRes[,"evenness"]))$p.value
          return( rbind(data.frame(cohort = x, tmp.FC),
                        data.frame(cohort=x, metric = "deltaEvenness", FC = deltaEvenness, Pvalue = Pvalue)))
}))

# CSR
csrFC = getFC(CSRData, metricName = "CSR")
shmFC = getFC(SHMData, metricName = "shm")
allIgProp$IgG1_plus_IgG3 = allIgProp$IGHG1 + allIgProp$IGHG3
IgG1FC = getFC(allIgProp, metricName = "IgG1_plus_IgG3")
# PDL1
PD1_FC = getFC( allTPMdata, metricName = "PDCD1")
PDL1_FC = getFC( allTPMdata, metricName = "CD274")
CTLA4_FC = getFC( allTPMdata, metricName = "CTLA4")
IGHC_FC = getFC( allTPMdata, metricName = "IGHC")

load("../htseq_vdj_from_David/immuneDeconv/All_FC_of_immuneDeconv_in_response.Rdata")
#Bcell = allFCdata[grepl("^B cell_|T cell CD8+", allFCdata$cellType),]
Bcell = allFCdata[grepl("MCP_COUNTER", allFCdata$cellType),]
Bcell.sub = Bcell[,c(2,1,3,4)]
colnames(Bcell.sub)[2] = "metric"

allFC = rbind(cpkFC, shmFC, IgG1FC, csrFC, Bcell.sub, PD1_FC, PDL1_FC, CTLA4_FC, IGHC_FC)
allFC$log2FC = log2(allFC$FC)
summary(allFC$log2FC)
allFC$log2FC[ allFC$log2FC > 2 ] = 2
allFC$log2FC[ allFC$log2FC < -2 ] = -2
#allFC = allFC %>% group_by(cohort) %>% mutate(normLog2FC = scale_rows(log2FC))
summary(allFC$normLog2FC)
allFC$Padj = p.adjust(allFC$Pvalue, method="BH")
allFC$cohort = sapply(allFC$cohort, function(x){ renameCohort(x, cohortNum = cohortName) })
unique(allFC$cohort)
# sort cohorts by Trust4 report reads
T4read = subset(allFC, metric == "normT4CDR3")
T4read = T4read[order(T4read$FC,decreasing = T),]
T4read$cohortname = name2cohort[ T4read$cohort ]
T4read.Pre = T4read$cohort[!grepl("On|EDT", T4read$cohortname)]
T4read.On = T4read$cohort[grepl("On|EDT", T4read$cohortname)]
cohortT4Order = c(as.character(T4read.Pre), as.character(T4read.On))
save(cohortT4Order, file = "cohort_order_by_Trust4_IGH_abundance.Rdata")

allFC$cohort = factor(allFC$cohort, levels =cohortT4Order)
allFC$metric = sapply(as.character(allFC$metric), function(x){
    if( x == "CPK"){
      return("Diversity")
    }else if(x=="normUniqCDR3"){
      return("Richness")
    }else if(x=="shm"){
      return("SHM")
    }else if(x=="IgG1_plus_IgG3"){
      return("(IgG1+IgG3)%")
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

#========fold change==========
subFC = subset(allFC, metric %in% c("(IgG1+IgG3)%","Richness","Diversity","TRUST4_IGHabundance","B cell_MCP_COUNTER","B cell_XCELL","T cell CD8+_MCP_COUNTER","PDL1","PD1","CTLA4") )#"deltaEvenness","CSR","SHM","IGHC"
subFC = subset(subFC, !grepl("SQ|LN|bladder|kidney|lymph|ureter",cohort) )
subFC = subset(subFC, !grepl("Gide_PD1|Gide_ipiPD1",cohort) )
subFC = subset(subFC, !grepl("Riaz\\(",cohort) )
subFC$Padj = p.adjust(subFC$Pvalue, method = "BH")
subFC.exDeltaEven = subset(subFC, metric != "deltaEvenness")
subFC.exDeltaEven$metric = factor(subFC.exDeltaEven$metric, levels = c("CSR","SHM","(IgG1+IgG3)%","Richness","Diversity","TRUST4_IGHabundance","B cell_MCP_COUNTER","B cell_XCELL","T cell CD8+_MCP_COUNTER","PDL1","PD1","CTLA4","IGHC"))
subFC.exDeltaEven$Padj = p.adjust(subFC.exDeltaEven$Pvalue, method = "BH")
pdf("Fig_all_metric_in-heatmap.pdf", height = 8, width = 10)
source("~/Dropbox (Partners HealthCare)/Immunotherapy/src/heatmap_size.R")
p1 <- corPlot(subFC.exDeltaEven, x_axis = "cohort", y_axis = "metric", titleName = "Responder / nonResponder",
        rho = "log2FC", pvalue = "Padj", fdr_cutoff = 0.2, theme = mytheme_noXlab,
        min_value = -2, max_value = 2, bar_anno = c(-2,0,2))
subFC.exDeltaEven = subFC.exDeltaEven[ !is.na(subFC.exDeltaEven$cohort), ]
corPlot(subFC.exDeltaEven, x_axis = "cohort", y_axis = "metric", titleName = "Responder / nonResponder",
        rho = "log2FC", pvalue = "Padj", fdr_cutoff = 0.1, theme = mytheme_grid,
        min_value = -2, max_value = 2, mean_value = 0,bar_anno = c(-2,0,2))
dev.off()

save( renameCohort, cohortT4Order, cohortName, my_theme_nocolor,file = "CohortT4Order_and_name_for_plot.Rdata")
