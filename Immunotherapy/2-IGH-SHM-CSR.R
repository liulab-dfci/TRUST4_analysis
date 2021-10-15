setwd("/Users/zyouyang/Dropbox (Partners HealthCare)/Immunotherapy/inter_cohorts")
library(dplyr)
library(ggpubr)
library(data.table)
source("/Users/zyouyang/Dropbox (Partners HealthCare)/Immunotherapy/src/heatmap_size.R")
load("cohort_color_for_plot.Rdata")
seqMisPos = function(seqs){
  if(length(seqs)==1){0
    return(NA)
  }else{
    test = matrix(unlist(strsplit(seqs, "")), nrow = length(seqs), byrow = T);
    count_n = apply(test, 2, function(x) table(x)[c("A", "T", "C", "G", "-")]);
    pos = which(colSums(sign(count_n[-5,]), na.rm=T)>1);
    return(paste(pos, collapse=","))
  }
}

#---------analysis----------
cohorts = list.files(path = ".",pattern = "igh\\.Rdata$")
cohorts
selectIg1 = clusterSizeDo = "N"
SHMdo = "Y"
VgeneDo = "N"
SHMData = NULL
CSRData = NULL
Ig1Switch = NULL
clusterSize = NULL
cdr3lenData = NULL
cdr3lenByIgData = NULL
Vusage = NULL
Jusage = NULL
for(i in cohorts){
  #i = cohorts[8]
  i = gsub(pattern = "_igh.Rdata",replacement = "",i)
  igh = get(load(paste0(i,"_igh.Rdata")))
  igh = subset(igh, is_complete=="Y" & !grepl("_",CDR3_aa))
  clinicalFile = list.files(path = ".",pattern = paste0(i,".*.","clinical.Rdata"))
  for(j in clinicalFile){
    #j = clinicalFile[1]
    print(j)
    clinical = get(load(j))
    clinical = clinical[!is.na(clinical$Response),]
    clinical$responseLabel= sapply(clinical$Response,function(x){if(x==1){return("responder")}else{return("non-responder")}})
    clinical$responseLabel = factor(clinical$responseLabel,levels = c("responder","non-responder"))
    # IGH SHM
    if( SHMdo == "Y"){
      # filter low freq CDR3
      rerun = "N"
      if( rerun == "Y"){
        igh = igh %>% group_by(Patient,clusterID) %>% mutate(maxCount = max(count))
        igh = subset(igh, count > maxCount*0.001)
        mismatch = igh[,c('Patient','clusterID','CDR3_DNA')] %>% group_by(Patient,clusterID) %>%filter(n_distinct(CDR3_DNA)>1) %>% summarise(size=length(CDR3_DNA), length=max(nchar(CDR3_DNA)), misPos=seqMisPos(CDR3_DNA))
        mismatch$misCount <- sapply(mismatch$misPos,function(x){if(is.na(x)){return(0)} else{ posList=strsplit(as.character(x),split=','); return(length(posList[[1]]))}})
        SHM <- mismatch %>% group_by(Patient) %>% summarise(shm = sum(misCount)/sum(length))
        #summary(SHM$shm)
        SHMinfo = merge(SHM,clinical,by="Patient")
        save(mismatch,SHMinfo, file = gsub("_clinical.Rdata","_SHM.filterLowFreqCDR3.Rdata",j))
      }

      load(gsub("_clinical.Rdata","_SHM.filterLowFreqCDR3.Rdata",j))
      SHMData = rbind(SHMData, data.frame(SHMinfo[,c("Patient","shm","responseLabel")],cohort=gsub("_clinical.Rdata","",j)))
    }
    # class switch
    igh_rmUnknownC = subset(igh, C!="*")
    uniqCDR3 = igh_rmUnknownC %>% group_by(Patient) %>% summarise(CDR3size = n_distinct(clusterID))
    igh_cluster = igh_rmUnknownC %>% group_by(Patient,clusterID) %>% filter(length(unique(C))>1) %>% summarise(newC = paste(unique(C),collapse=","),uniqCnum = length(unique(C)))
    # cluster size
    if( clusterSizeDo == "Y"){
      igh_cluster_size = igh_rmUnknownC %>% group_by(Patient,clusterID) %>% summarise(clusterSize = n())
      igh_cluster_size = merge(igh_cluster_size, clinical[,c("Patient","responseLabel")])
      clusterSize = rbind(clusterSize, data.frame(igh_cluster_size,cohort=gsub("_clinical.Rdata","",j)))
    }
   ## CSR = number_of_clusters_with_CS/number_of_unique_CDR3
    clusterNum = igh_cluster %>% group_by(Patient) %>% summarise(clusterNum = n_distinct(clusterID))
    CSR = Reduce(function(df1,df2) merge(df1,df2,by="Patient",all.x=T), list(clinical[,c("Patient","responseLabel")],uniqCDR3, clusterNum))
    if( nrow(CSR)==0 ){ next }
    CSRData = rbind(CSRData, data.frame(CSR,cohort=gsub("_clinical.Rdata","",j)))
    # IGH1 switch level
    if( selectIg1 == "Y"){
      igh1cs = igh_cluster[grep('IGHG1',igh_cluster$newC),]
      igh_rmUnknownC$label = sapply(igh_rmUnknownC$clusterID, function(x) { if(x %in% igh_cluster$clusterID){return("clustered")}else{return("not-clustered")}})
      igh_rmUnknownC$ig1.label = sapply(igh_rmUnknownC$clusterID, function(x) { if(x %in% igh1cs$clusterID){return("clustered")}else{return("not-clustered")}})
      igh1switchProp = igh_rmUnknownC %>% group_by(Patient,ig1.label) %>% summarise(ig1switch = n())
      igh1switchProp = igh1switchProp %>% group_by(Patient) %>% mutate(freq=ig1switch/sum(ig1switch))
      igh1switchProp = merge(igh1switchProp, clinical[,c("Patient","responseLabel")],by="Patient")
      Ig1Switch = rbind(Ig1Switch,data.frame(igh1switchProp, cohort=gsub("_clinical.Rdata","",j)))
    }
    # CDR3 length
    cdr3len = igh %>% group_by(Patient) %>% summarise(medianLen = median(CDR3aalen))
    cdr3len = merge(cdr3len, clinical[,c("Patient","responseLabel")],by="Patient")
    cdr3lenData = rbind(cdr3lenData, data.frame(cdr3len,cohort=gsub("_clinical.Rdata","",j)))
    # CDR3 length by isotype
    cdr3len = igh %>% group_by(Patient,C) %>% summarise(medianLen = median(CDR3aalen))
    cdr3len = merge(cdr3len, clinical[,c("Patient","responseLabel")],by="Patient")
    cdr3lenByIgData = rbind(cdr3lenByIgData, data.frame(cdr3len,cohort=gsub("_clinical.Rdata","",j)))
    if( VgeneDo == "Y"){
        igh = igh %>% group_by(Patient) %>% mutate(freq = count/sum(count))
        # V gene
        igh$Vgene = sapply(strsplit(igh$V,"\\*"),"[",1)
        vgene = igh %>% group_by(Patient, Vgene) %>% summarise(freq=sum(freq))
        vgene = merge(vgene, clinical[,c("Patient","responseLabel")],by="Patient")
        Vusage = rbind(Vusage,data.frame(vgene,cohort=gsub("_clinical.Rdata","",j)))
        # J gene
        igh$Jgene = sapply(strsplit(igh$J,"\\*"),"[",1)
        jgene = igh %>% group_by(Patient, Jgene) %>% summarise(freq=sum(freq))
        jgene = merge(jgene, clinical[,c("Patient","responseLabel")],by="Patient")
        Jusage = rbind(Jusage,data.frame(jgene,cohort=gsub("_clinical.Rdata","",j)))
    }
  }
}

CSRData = CSRData[!is.na(CSRData$responseLabel) & !is.na(CSRData$CDR3size),]
CSRData[is.na(CSRData$clusterNum),"clusterNum"] = 0
CSRData$CSR = CSRData$clusterNum/CSRData$CDR3size
#CSRData$cohort = factor(CSRData$cohort,levels = cohortOrder)

save(SHMData,file="SHM_of_cohorts.IGH.filterLowFreqCDR3.Rdata")
save(CSRData,file="CSR_of_cohorts.IGH.Rdata")
save(cdr3lenData, file = "CDR3_AA_length.IGH.Rdata")
save(Vusage,Jusage, file = "VJgene_usage.Rdata")

##---------plot------------
load("SHM_of_cohorts.IGH.filterLowFreqCDR3.Rdata")
load("CSR_of_cohorts.IGH.Rdata")
load("CDR3_AA_length.IGH.Rdata")
load("cohort_color_for_plot.Rdata")
load("CohortT4Order_and_name_for_plot.Rdata")

# SHM
SHMData$cohort = sapply(SHMData$cohort, function(x)( return( renameCohort(x, cohortName))))
SHMData$cohort = factor(SHMData$cohort, levels=cohortT4Order)
SHM.Subdata = subset(SHMData, shm<=0.2)
SHM.Subdata = subset(SHM.Subdata, !grepl("LN|SQ|kidney|bladder|lymph|ureter", cohort) & !grepl("Gide_PD1|Gide_ipiPD1", cohort) & !grepl("Riaz\\(", cohort))
pdf("Fig_SHM_of_IGH.pdf")
ggboxplot(SHM.Subdata, x = "cohort", y = "shm", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),ylim=c(0,0.2), ylab = "SHM ratio") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
dev.off()


# CSR
CSRData$cohort = sapply(CSRData$cohort, function(x)( return( renameCohort(x, cohortName))))
CSRData$cohort = factor(CSRData$cohort, levels=cohortT4Order)
CSRData = subset(CSRData, !grepl("LN|SQ|kidney|bladder|lymph|ureter", cohort) & !grepl("Gide_PD1|Gide_ipiPD1", cohort) & !grepl("Riaz\\(", cohort))
CSRData.sub = subset(CSRData, CDR3size>0)
pdf("Fig_CSR_of_IGH.pdf")
ggboxplot(CSRData.sub, x = "cohort", y = "CSR", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),title = "CDR3size>=5", ylab = "cluster switch ratio") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# scale y max 0.5
CSRData.sub = subset(CSRData, CSR<0.5)
ggboxplot(CSRData.sub, x = "cohort", y = "CSR", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),ylim=c(0,0.2)) +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
dev.off()
