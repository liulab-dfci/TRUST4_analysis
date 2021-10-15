setwd('~/Dropbox (Partners HealthCare)/Immunotherapy/inter_cohorts/')
library(ggpubr)
library(dplyr)
library(gridExtra)
Ig.cols=colors()[c(26,652,3,47,139,419,116,120,115,23)] #47 -> 68
Ig.genes=c('IGHD','IGHM','IGHE','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4')
names(Ig.cols)=c(Ig.genes,'Unidentified')

#-------function ---------
entropy <- function(data, norm=F, cutoff = 0.95){
  if(length(data)==0){
    return(NA)
  }
  #data = c(1:10)
  #data = pt44$count
  data = sort( prop.table(data), decreasing = T)
  #print( length(data) )
  res = -sum(data*log2(data))
  if(norm){
    # top 95% abundant clones for evenness calculation
    tmp.cumsum = cumsum(data)
    nStart = min( which(tmp.cumsum >= cutoff))
    data = data[1:nStart]
    res = -sum(data*log2(data))
    return(res/log2(length(data)))
  }
  else{
    return(res)
  }
}

entropyFilterLowAbudant <- function(data, norm=F, freq= 0.01){
  # filter out clones with freq < 0.01
  data = prop.table( data )
  data = data[ data >= freq ]
  data = prop.table( data )
  res = -sum(data*log2(data))
  if(norm){
    return(res/log2(length(data)))
  }
  else{
    return(res)
  }
}
topFreq <- function(data,percentage=0.02){
  # sum frequency of top 1% clusters
  if(length(data)< 1/percentage){
    return(NA)
  }
  data = sort(data,decreasing = T)
  n = as.integer( length(data)*percentage )
  res = sum(data[1:n])/sum(data)
  return(res)
}

largeClonesN <- function(data, cutoff = 0.005){
  data = prop.table(data)
  nLarge = length( which( data >= cutoff) )
  return( nLarge )
}
largeClonesProp <- function(data, cutoff = 0.005){
  data = prop.table(data)
  nLarge = length( which( data >= cutoff) )
  return( nLarge/length(data) )
}
#---------analysis----------
cohorts = list.files(path = ".",pattern = "igh\\.Rdata$")
cohorts
allData = NULL
allIgCount = NULL
allIgProp = NULL
allTop5Freq = NULL
allWtMean = NULL
clusterSizeAsAbudant = "N"
for(i in cohorts){
  #i = cohorts[9]
  i = gsub(pattern = "_igh.Rdata",replacement = "",i)
  clinicalFile = list.files(path = ".",pattern = paste0(i,".*.","clinical.Rdata"))
  #igh.partial = get(load(paste0(i,"_igh.partial.Rdata")))
  #igh.partial = subset(igh, !grepl("_",CDR3_aa))
  igh = get(load(paste0(i,"_igh.Rdata")))
  igh = subset(igh, !grepl("_",CDR3_aa) )
  #igh = subset(igh, count>1)
  seqDepth = get(load(paste0(i,"_seqDepth.Rdata")))
  #vdjCov = get(load(paste0(i,"_vdjCov.Rdata")))
  for(j in clinicalFile){
    #j = clinicalFile[1]
    print(j)
    clinical = get(load(j))
    clinical = clinical[!is.na(clinical$Response),]
    clinical$responseLabel= sapply(clinical$Response,function(x){if(x==1){return("responder")}else{return("non-responder")}})
    # benefit or not
    #if( grepl("Riaz", i)){
    #  clinical$responseLabel= sapply(clinical$RECIST,function(x){if(x==0){return("non-responder")}else{return("responder")}})
    #}
    clinical$responseLabel = factor(clinical$responseLabel,levels = c("responder","non-responder"))
    # complete CDR3
    igh.complete = subset(igh, is_complete == 'Y' & CDR3_aa !="out_of_frame")
    #============= wighted mean similarity ================
    ighSim = igh.complete %>% group_by(Patient) %>% summarise(wtMean = weighted.mean(similarity, count, na.rm=T))
    ighSimInfo = merge(ighSim, clinical[,c("Patient","responseLabel")], by = "Patient")
    allWtMean = rbind(allWtMean, data.frame(ighSimInfo,cohort=gsub("_clinical.Rdata","",j)))
    #================ CPK ================
    ighCluster = igh.complete %>% group_by(Patient,clusterID) %>% summarise(count=sum(count),clusterSize=n_distinct(CDR3_DNA))
    ighCluster = subset(ighCluster, count>1)
    if( clusterSizeAsAbudant == "Y"){
      ighCluster = ighCluster %>% mutate(count = clusterSize)
    }
    #ighStat = ighCluster %>% group_by(Patient) %>%  summarise(libsize = sum(count),uniqCDR3 = n_distinct(clusterID),entropy = entropy(count,cutoff = 0.9999), evenness = entropy(count,norm=T),t5freq = topFreq(count, percentage = 0.02) ) %>% mutate(CPK = (uniqCDR3/libsize)*1000, clonality = 1-evenness)
    ighStat = ighCluster %>% group_by(Patient) %>%  summarise(minClusterSize = min(clusterSize), maxClusterSize = max(clusterSize), libsize = sum(count),uniqCDR3 = n_distinct(clusterID),entropy = entropy(count), evenness = entropy(count,norm=T,cutoff = 0.95),t2freq = topFreq(count, percentage = 0.02), nLarge = largeClonesN(count), propLarge = largeClonesProp(count) ) %>% mutate(CPK = (uniqCDR3/libsize)*1000, clonality = 1-evenness)
    #ighStat = na.omit(ighStat)
    #ighLibsize = igh.partial %>% group_by(Patient) %>% summarise(t4libsize = sum(count))
    #ighStat = ighStat %>% left_join(ighLibsize,by="Patient") %>% mutate(CPK=(uniqCDR3/t4libsize)*1000)
    #ighInfo = Reduce(function(df1,df2) merge(df1,df2,by="Patient"),list(ighStat,clinical[,c("Patient","responseLabel")],vdjCov)) %>% mutate(CPK = (uniqCDR3/coverage)*1000)
    ighInfo = Reduce(function(df1,df2) merge(df1,df2,by="Patient"),list(ighStat,clinical[,c("Patient","responseLabel")],seqDepth)) #%>% mutate(CPK = (uniqCDR3/totalReads)*1000)
    allData = rbind(allData, data.frame(ighInfo[,c(colnames(ighStat),"totalReads","responseLabel")],cohort=gsub("_clinical.Rdata","",j)))

    # Ig abundance
    for (Ig in Ig.genes) {
      #print(Ig)
      IgG1 = igh %>% group_by(Patient) %>% filter(grepl(Ig,C)) %>% summarise(Cnum=sum(count))
      IgInfo = merge(IgG1, clinical[,c('Patient','responseLabel')], by="Patient")
      IgNormalize = merge(IgInfo,seqDepth[,c('Patient','totalReads')],by="Patient")
      if(nrow(IgNormalize)>0){
        allIgCount = rbind(allIgCount, data.frame(IgNormalize,Cgene=Ig,cohort=gsub("_clinical.Rdata","",j)))
      }
    }
    # Ig fraction
    all.tt = table(igh$Patient)
    Ig.mat=matrix(0,ncol=length(all.tt),nrow=10,dimnames = list(names(Ig.cols),
                                                                names(all.tt)))
    for(Ig in names(Ig.cols)){
      ## Extract samples with the Ig isotype
      tmp.vv=grep(Ig,igh[,'C'])
      ## Count the amount of CDR3 with the Ig isotype for previous samples
      tmp.tt=table(igh[tmp.vv,'Patient'])
      ## Get the ratio of the Ig isotyoe by dividing the total amount of CDR3
      Ig.mat[Ig,names(tmp.tt)]=tmp.tt
    }

    Ig.mat = as.data.frame(t(Ig.mat))
    Ig.mat$total = apply(Ig.mat,1,sum)
    Ig.mat = subset(Ig.mat, total>0)

    ##Ig proportion
    IgProp = prop.table(as.matrix(Ig.mat[,1:9]),margin=1)
    IgProp = as.data.frame(na.omit(IgProp))
    IgPropInfo = merge(IgProp, clinical[,c('Patient','responseLabel')], by.x=0,by.y="Patient")
    allIgProp = rbind(allIgProp,data.frame(IgPropInfo,cohort=gsub("_clinical.Rdata","",j)))
  }
}
save(allData,allIgProp,allIgCount, file="CPK_allData.Rdata")

#============use Trust4 IGH abundance order======
load("cohort_order_by_Trust4_IGH_abundance.Rdata")
cohortOrder = cohortT4Order

my_theme = theme(legend.title = element_text(size = 7),
                 legend.text = element_text(size=5),
                 legend.key.width=unit(1, "lines"),
                 plot.title = element_text(size=10, face="bold", hjust=0.5),
                 plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),
                 axis.text.x = element_text(angle = 90, size=9,face="bold",
                                            #color = c(rep("deepskyblue",length(effectCohorts)), rep("black",length(nonEffectCohorts)-1),"red3")),
                                            #color = cohortLabelColor[cohortOrder]),
                                            color = cohortLabelColor[cohortT4Order]),
                 axis.title.x=element_blank(),
                 axis.text.y = element_text(face="bold"),
                 panel.background = element_rect(fill = NA, colour = "black"),
                 panel.grid.major = element_blank())# linetype="dashed"

my_theme_nocolor = theme(legend.title = element_text(size = 7),
                 legend.text = element_text(size=5),
                 legend.key.width=unit(1, "lines"),
                 plot.title = element_text(size=10, face="bold", hjust=0.5),
                 plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),
                 axis.text.x = element_text(angle = 90, size=9,face="bold",
                                            color = "black"),
                 axis.title.x=element_blank(),
                 axis.text.y = element_text(face="bold"),
                 panel.background = element_rect(fill = NA, colour = "black"),
                 panel.grid.major = element_blank())# linetype="dashed"
load("cohort_color_for_plot.Rdata")
cohortOrder = cohortT4Order
save(cohortOrder, cohortLabelColor, my_theme, my_theme_nocolor, file= "cohort_color_for_plot.Rdata")

#-----------plot------------------
load("CPK_allData.Rdata")
table(allData$cohort, allData$responseLabel)
allData$cohort = sapply(allData$cohort, function(x) { return( renameCohort(x, cohortName))})
allData$cohort = factor(allData$cohort, levels=cohortT4Order)
allData$time = sapply(allData$cohort, function(x) {
    if(grepl("On|EDT",x)){
      return("On")
    }else{
      return("Pre")
    }
})
allData$time = factor(allData$time, levels=c("Pre","On"))
summary(allData$uniqCDR3)
cutoff= 0
allData.sub = subset(allData, uniqCDR3 > cutoff)
allData.sub = subset(allData.sub, !grepl("LN|SQ|kidney|bladder|lymph|ureter", cohort))
allData.sub = subset(allData.sub, !grepl("Gide_PD1|Gide_ipiPD1", cohort))
allData.sub = subset(allData.sub, !grepl("Riaz\\(", cohort))
pdf(paste0("Fig_igh_entropy_response.complete.oneLine.cluster.uniqCluster",cutoff,"XCELL.pdf"))
# entropy
ggboxplot(allData.sub, x = "cohort", y = "entropy", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"))+
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# unique CDR3
ggboxplot(allData.sub, x = "cohort", y = "uniqCDR3", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold")) +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# normalized unique CDR3
allData.sub$normUniqCDR3 = (allData.sub$uniqCDR3/allData.sub$totalReads)*1000000
allData.sub.exGC = subset(allData.sub, !grepl("STAD",cohort) )
ggboxplot(allData.sub.exGC, x = "cohort", y = "normUniqCDR3", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),ylab="normUniqCluster") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# mGC
allData.sub.Kim = subset(allData.sub, grepl("STAD",cohort))
ggboxplot(allData.sub.Kim, x = "cohort", y = "normUniqCDR3", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),ylab="normUniqCluster") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
#grid.arrange(p1,p2,nrow=1)
# evenness
ggboxplot(allData.sub, x = "cohort", y = "evenness", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold")) +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# CPK
# non-responder outlier
ggboxplot(allData.sub, x = "cohort", y = "CPK", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),title="complete CPK") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# exclude Kim and Hugo
ggboxplot(allData.sub.exGC, x = "cohort", y = "CPK", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),title="complete CPK") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# Kim and Hugo
ggboxplot(allData.sub.KimHugo, x = "cohort", y = "CPK", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),title="complete CPK") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor

allData.sub$CPKn = (allData.sub$uniqCDR3/allData.sub$libsize)*1000
ggboxplot(allData.sub, x = "cohort", y = "CPKn", colsor = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),title="completeCPK") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# top 2% frequency
ggboxplot(allData.sub, x = "cohort", y = "t2freq", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),ylab = "frequency of top 2% clonotypes") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# weighted mean similarity
allWtMean$cohort = factor(allWtMean$cohort, levels = cohortOrder)
ggboxplot( allWtMean, x = "cohort", y = "wtMean", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
           font.label = list(size = 9, face = "bold"),ylab = "weighted mean similarity") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
dev.off()

# Ig normalized abundance
pdf("Fig_normalized_Ig_abundance.cluster.pdf")
allIgCount$normCnum = (allIgCount$Cnum/allIgCount$totalReads)*1000000
for(Ig in Ig.genes){
  subdata = subset(allIgCount, Cgene==Ig)
  print( ggboxplot(subdata, x = "cohort", y = "normCnum", color = "responseLabel",add="jitter",add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
            font.label = list(size = 9, face = "bold"),title = Ig, ylab = "normalized count") +
    stat_compare_means(label =  "p.signif",aes(group=responseLabel))+my_theme )
}
dev.off()

allIgProp$cohort = sapply(allIgProp$cohort, function(x) { return( renameCohort(x, cohortName))})
allIgProp$cohort = factor(allIgProp$cohort, levels=cohortT4Order)
allIgProp= subset(allIgProp, !grepl("LN|SQ|kidney|bladder|lymph|ureter", cohort))
allIgProp = subset(allIgProp, !grepl("Gide_PD1|Gide_ipiPD1", cohort))
allIgProp = subset(allIgProp, !grepl("Riaz\\(", cohort))
pdf("Fig_Ig_fraction.oneLine.cluster.xcell.pdf")
for(Ig in Ig.genes){
  print( ggboxplot(allIgProp, x = "cohort", y = Ig, color = "responseLabel",add="jitter",add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
                   font.label = list(size = 9, face = "bold"),title = Ig) +
           stat_compare_means(label =  "p.signif",aes(group=responseLabel))+my_theme_nocolor )
}
allIgProp$IgA1divIgG1 = allIgProp$IGHA1/allIgProp$IGHG1
ggboxplot(allIgProp, x = "cohort", y = "IgA1divIgG1", color = "responseLabel",add="jitter",add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),title = "IgA1_div_IgG1", ylab = "IgA1_div_IgG1") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel))+my_theme_nocolor
allIgProp$IgG1_plus_G3 = allIgProp$IGHG1+allIgProp$IGHG3
ggboxplot(allIgProp, x = "cohort", y = "IgG1_plus_G3", color = "responseLabel",add="jitter",add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),title = "IgG1_plus_G3", ylab = "(IgG1+IgG3)%") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel))+my_theme_nocolor
allIgProp$IgA1_plus_A2 = allIgProp$IGHA1+allIgProp$IGHA2
ggboxplot(allIgProp, x = "cohort", y = "IgA1_plus_A2", color = "responseLabel",add="jitter",add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold"),title = "IgA1_plus_A2", ylab = "IgA1_plus_A2") +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel))+my_theme_nocolor
dev.off()

#=============trb===============
allData = NULL
for(i in cohorts){
  #i = cohorts[8]
  i = gsub(pattern = "_igh.Rdata",replacement = "",i)
  clinicalFile = list.files(path = ".",pattern = paste0(i,".*.","clinical.Rdata"))
  trb = get(load(paste0(i,"_trb.Rdata")))
  for(j in clinicalFile){
    #j = clinicalFile[1]
    print(j)
    clinical = get(load(j))
    clinical = clinical[!is.na(clinical$Response),]
    clinical$responseLabel= sapply(clinical$Response,function(x){if(x==1){return("responder")}else{return("non-responder")}})
    clinical$responseLabel = factor(clinical$responseLabel,levels = c("responder","non-responder"))
    # complete CDR3
    trb.complete = subset(trb, is_complete == 'Y' & CDR3_aa !="out_of_frame")
    trbCluster = trb.complete %>% group_by(Patient, clusterID) %>% summarise(count=sum(count),freq=sum(freq))
    trbStat = trbCluster %>% group_by(Patient) %>%  summarise(libsize = sum(count),uniqCDR3 = n_distinct(clusterID),entropy = entropy(count), evenness = entropy(count,norm=T)) %>% mutate(CPK = (uniqCDR3/libsize)*1000, clonality = 1-evenness)
    trbStat = na.omit(trbStat)
    trbInfo = merge(trbStat,clinical,by="Patient")
    allData = rbind(allData, data.frame(trbInfo[,c(colnames(trbStat),"responseLabel")],cohort=gsub("_clinical.Rdata","",j)))
  }
}
save(allData, file = "TCR_metric.Rdata")

pdf("Fig_trb_entropy_response.complete.oneLine.pdf")
allData$cohort = factor(allData$cohort, levels=cohortOrder)
# entropy
ggboxplot(allData, x = "cohort", y = "entropy", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold")) +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
# evenness
ggboxplot(allData, x = "cohort", y = "evenness", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold")) +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
ggplot(allData, aes(x = cohort, y = evenness, color = responseLabel, fill = responseLabel)) +
  geom_split_violin() + my_theme_nocolor

# CPK
ggboxplot(allData, x = "cohort", y = "CPK", color = "responseLabel",add="jitter",width=0.5,add.params = list(size = 0.5, jitter = 0.2), xlab="",font.title = list(size=10),
          font.label = list(size = 9, face = "bold")) +
  stat_compare_means(label =  "p.signif",aes(group=responseLabel)) + my_theme_nocolor
dev.off()
