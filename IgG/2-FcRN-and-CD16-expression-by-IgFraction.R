library(ggpubr)
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggsci)
library( ppcor )
my_theme = theme(legend.title = element_blank(),
                 legend.text = element_text(size=17, face = "bold"),
                 legend.key.width=unit(1, "lines"),
                 plot.title = element_text(size=10, face="bold", hjust=0.5),
                 #plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),
                 axis.text.x = element_text(angle = 90, size=17,face="bold"),
                 axis.title.x=element_blank(),
                 axis.text.y = element_text(size = 17, face="bold"),axis.title.y = element_text(size = 17,face="bold"),
								 #legend.position="right",
                 panel.grid.major = element_blank())
#setwd("/liulab/zouyang/TCGA/IGH/1_SHM")
setwd("/liulab/lsong/projects/TRUST4_Ouyang_analysis/Fig5")
source('/homes/zouyang/Immunotherapy/src/heatmap_size.R')
source("/liulab/lsong/projects/TRUST4_Ouyang_analysis/func_survival.R")
source("/liulab/lsong/projects/TRUST4_Ouyang_analysis/func_survival_4class.R")

#================combine FcRN expression with IgG1 or IgA1 CPK==============
load("/liulab/zouyang/TCGA/IGH/1_SHM/sampleInfo_and_cancer_type.Rdata")
tumorPurity=get(load('/liulab/lsong/data/bcr/BCR_NG_TRUST/data/TRUST3/tumorPurity.Rdata'))
exprDat = get(load('/liulab/zouyang/TCGA/QuantFile/TPM_of_TCGA.proteinCoding.Rdata'))
idIndex = colnames(exprDat)
#-----C1Q expression------------
C1Q = apply(exprDat[c("C1QA", "C1QB", "C1QC"),],2,mean)

##=========IgG1 + igG3 fraction==========
load("/liulab/zouyang/TCGA/IGH/1_SHM/Ig_mat.Rdata")
Ig.matProp = prop.table(as.matrix(Ig.mat[,3:11]),margin=1)
Ig.matProp = as.data.frame(Ig.matProp)
rownames(Ig.matProp) = Ig.mat$TCGA_id
Ig.mat = Ig.matProp
Ig.mat$IgG1_plus_IgG3 = Ig.mat[,'IGHG1']+Ig.mat[,'IGHG3']
Ig.mat$IGHA = Ig.mat[,'IGHA1']+Ig.mat[,'IGHA2']
Ig.mat$IGHG = Ig.mat[,'IGHG1'] + Ig.mat[,'IGHG2'] + Ig.mat[,'IGHG3'] + Ig.mat[,'IGHG4']

# IgG1 + IgG3
igg1 = Ig.mat$IgG1_plus_IgG3
names(igg1) = rownames(Ig.matProp)

combineFcRN = data.frame(expr = unlist(exprDat['FCGRT', idIndex]), igg1count = igg1[ idIndex ],
              purity = tumorPurity[ substr(idIndex, 1 , 16)], cancer = cancer.types[ idIndex ],
              CD16a = unlist(exprDat['FCGR3A', idIndex]), C1Q = C1Q[ idIndex ])

save( combineFcRN, file = "combine_data_of_FcRN_and_IgG1_fraction_and_C1Q.Rdata")

my_comparisons <- list(c("low FcRN + low (IgG1+IgG3)%", "low FcRN + high (IgG1+IgG3)%"), c("high FcRN + low (IgG1+IgG3)%", "high FcRN + high (IgG1+IgG3)%"),
                      c("low FcRN + high (IgG1+IgG3)%", "high FcRN + high (IgG1+IgG3)%"),c("low FcRN + low (IgG1+IgG3)%", "high FcRN + low (IgG1+IgG3)%"))

#==============CD16a expression classified by FcRN and IgG1 fraction=============
fcrnIgG1group = Reduce( rbind, lapply( sort(unique(combineFcRN$cancer)), function(x){
      subdata = subset(combineFcRN, cancer == x)
      subdata = get2ClassBy37(subdata, first = "expr", second = "igg1count", cutoff=0.5)
      return(subdata)
}))
fcrnIgG1group$classLab = sapply(fcrnIgG1group$class, function(x){
       if(x == 1){
         return("low FcRN + low (IgG1+IgG3)%")
       }else if (x == 2){
         return("low FcRN + high (IgG1+IgG3)%")
       }else if (x == 3){
         return("high FcRN + low (IgG1+IgG3)%")
       }else if (x == 4){
         return("high FcRN + high (IgG1+IgG3)%")
       }})
fcrnIgG1group$FcrnLab = sapply(strsplit( as.character(fcrnIgG1group$classLab), " \\+ "), "[", 1)
fcrnIgG1group$FcrnLab = factor( fcrnIgG1group$FcrnLab, levels = c("low FcRN", "high FcRN"))
fcrnIgG1group$IgG1CPKLab = sapply(strsplit( as.character(fcrnIgG1group$classLab), " \\+ "), "[", 2)
fcrnIgG1group$IgG1CPKLab = factor(fcrnIgG1group$IgG1CPKLab,levels = c("low (IgG1+IgG3)%", "high (IgG1+IgG3)%"))
pdf("Fig_CD16a_in_FcRN_and_IgG1fraction.pdf",height = 17, width = 25)
fcrnIgG1group$log2CD16a = log2(fcrnIgG1group$CD16a)
fcrnIgG1group$log2C1Q = log2(fcrnIgG1group$C1Q)
fcrnIgG1group$classLab = factor(fcrnIgG1group$classLab, levels = c("low FcRN + low (IgG1+IgG3)%", "low FcRN + high (IgG1+IgG3)%", "high FcRN + low (IgG1+IgG3)%", "high FcRN + high (IgG1+IgG3)%"))
fcrnIgG1group$FcrnLab = factor(fcrnIgG1group$FcrnLab, levels = c("low FcRN", "high FcRN"))
fcrnIgG1group$IgG1CPKLab = factor(fcrnIgG1group$IgG1CPKLab, levels = c("low (IgG1+IgG3)%", "high (IgG1+IgG3)%"))
ggboxplot(fcrnIgG1group, x = "classLab", y = "log2CD16a", color = "IgG1CPKLab", ylab = expression(log[2]("CD16a")), facet.by = "cancer", nrow = 4, scales = "free_y") + stat_compare_means(label =  "p.signif", comparisons = my_comparisons) + my_theme
ggboxplot(fcrnIgG1group, x = "FcrnLab", y = "log2CD16a", color = "IgG1CPKLab", ylab = expression(log[2]("CD16a")), facet.by = "cancer", nrow = 4, scales = "free_y") + stat_compare_means(label =  "p.signif", aes(group = IgG1CPKLab)) + my_theme + scale_color_npg()
ggboxplot(fcrnIgG1group, x = "FcrnLab", y = "log2C1Q", color = "IgG1CPKLab", ylab = expression(log[2]("C1Q")), facet.by = "cancer", nrow = 4, scales = "free_y") + stat_compare_means(label =  "p.signif", aes(group = IgG1CPKLab)) + my_theme + scale_color_npg()
dev.off()


##========== NK cell by FcRN and IgG1CPK=============
fcrnIgG1group$barcode = substr( rownames(fcrnIgG1group), 1, 15)
load("/liulab/zouyang/TCGA/IGH/3_survival/immuneDeconv_result_of_TCGA.Rdata")
NKcell = immuneDecon[,grepl("NK cell|M1|M2|Macrophage", colnames(immuneDecon))]
#NKcell = immuneDecon[,grepl("MCP_COUNTER|XCELL", colnames(immuneDecon))]
NKcell = NKcell[, !grepl("CIBERSORT$", colnames(NKcell))]
FcRN_NKcell = merge(fcrnIgG1group, NKcell, by.x = "barcode", by.y = 0)
colnames(FcRN_NKcell)[14:ncol(FcRN_NKcell)] = gsub(" |-|/","_", colnames(FcRN_NKcell)[14:ncol(FcRN_NKcell)])

pdf("Fig_NKcell_in_FcRN_and_IgG1fraction.pdf",height = 20, width = 25)
FcRN_NKcell$log2MCP = log2(FcRN_NKcell$NK_cell_MCP_COUNTER)
ggboxplot(FcRN_NKcell, x = "FcrnLab", y = "log2MCP", ylab = expression(log[2]("NK_cell")),color = "IgG1CPKLab", facet.by = "cancer", nrow = 4,scales = "free_y") + stat_compare_means(aes(group = IgG1CPKLab),label =  "p.signif") + my_theme + scale_color_npg()
FcRN_NKcell$log2Macrophage = log2(FcRN_NKcell$Macrophage_Monocyte_MCP_COUNTER)
ggboxplot(FcRN_NKcell, x = "FcrnLab", y = "log2Macrophage", ylab = expression(log[2]("Macrophage_cell")),color = "IgG1CPKLab", facet.by = "cancer", nrow = 4,scales = "free_y") + stat_compare_means(aes(group = IgG1CPKLab),label =  "p.signif") + my_theme + scale_color_npg()
dev.off()


##===============tumor purity correlated with FcRN, CD16a, NK cell, macrophage cell infiltratoin===================
FcRN_NKcell_rmPurityNA = FcRN_NKcell[ !is.na(FcRN_NKcell$purity), ]
load("/liulab/zouyang/TCGA/IGH/1_SHM/metrics_of_TRB_and_IGH.Rdata")
igh_trb$barcode = substr( igh_trb$ID, 1, 15)
FcRN_NKcell_rmPurityNA = merge( FcRN_NKcell_rmPurityNA, igh_trb, by = "barcode" )
others = c("expr", "igg1count", "CD16a", colnames(igh_trb)[c(2:6,8:12)] ,colnames(FcRN_NKcell)[12:ncol(FcRN_NKcell)])
#others = c(colnames(FcRN_NKcell)[12:ncol(FcRN_NKcell)])
purityCorOther = Reduce( rbind, lapply( sort(unique(FcRN_NKcell_rmPurityNA$cancer)), function(x){
                    subdata = subset( FcRN_NKcell_rmPurityNA, cancer == x)
                    tmpDat = Reduce( rbind, lapply( others, function(y){
                              #print(x)
                              #print(y)
                              tmpCor = cor.test( subdata[, "purity"], subdata[,y], method = "spearman" )
                              if( y == "expr" ) {
                                  y = "FcRN"
                              }else if( y == "igg1count"){
                                  y = "(IgG1+IgG3)%"
                              }
                              return( data.frame( metric = y, corVal = tmpCor$estimate, pVal = tmpCor$p.value ) )
                    }))
                    return( data.frame( tmpDat, Disease = x) )
}))
purityCorOther$padj = p.adjust( purityCorOther$pVal, method = "BH" )
minVal = min(purityCorOther$corVal)
maxVal = max(purityCorOther$corVal)
pdf("Fig_cor_between_tumor_purity_and_FcRN_and_NKcell.pdf", width = 14, height = 14)
corPlot( purityCorOther, x_axis = 'Disease',y_axis = 'metric',rho = 'corVal',pvalue = 'padj',box_size = 5,
       max_value=maxVal, min_value=minVal, titleName = "cor between tumor purity and others ",
       bar_anno  = c( minVal, 0, maxVal) )
dev.off()

##===============expression of CD16a in all cancers=================
combineClass = get2ClassBy37(combineFcRN, first = "expr", second = "igg1count", cutoff=0.5)
combineClass$log2CD16a = log2(combineClass$CD16a)
combineClass$log2C1Q = log2(combineClass$C1Q)
combineClass$classLab = sapply(combineClass$class, function(x){
       if(x == 1){
         return("low FcRN + low (IgG1+IgG3)%")
       }else if (x == 2){
         return("low FcRN + high (IgG1+IgG3)%")
       }else if (x == 3){
         return("high FcRN + low (IgG1+IgG3)%")
       }else if (x == 4){
         return("high FcRN + high (IgG1+IgG3)%")
       }})
combineClass$classLab = factor(combineClass$classLab, levels = c("low FcRN + low (IgG1+IgG3)%", "low FcRN + high (IgG1+IgG3)%", "high FcRN + low (IgG1+IgG3)%", "high FcRN + high (IgG1+IgG3)%"))
combineClass$FcrnLab = sapply(strsplit( as.character(combineClass$classLab), " \\+ "), "[", 1)
combineClass$FcrnLab = factor(combineClass$FcrnLab,levels = c("low FcRN", "high FcRN"))
combineClass$IgG1CPKLab = sapply(strsplit( as.character(combineClass$classLab), " \\+ "), "[", 2)
combineClass$IgG1CPKLab = factor(combineClass$IgG1CPKLab,levels = c("low (IgG1+IgG3)%", "high (IgG1+IgG3)%"))
pdf("Fig_CD16a_in_FcRN_and_IgG1fraction.allCancer.pdf", width = 3, height = 10)
mytheme = theme(legend.title = element_blank(),
           legend.text = element_text(size=17, face = "bold"),
           legend.key.width=unit(1, "lines"),
           plot.title = element_text(size=17, face="bold", hjust=0.5),
           #plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),
					 legend.direction='vertical',
           axis.text.x = element_text(angle = 90, size=17,face="bold"),
           axis.title.x=element_blank(),
           axis.text.y = element_text(size=17, face="bold"),axis.title.y = element_text(size=17,face="bold"),
           panel.grid.major = element_blank())
ggboxplot(combineClass, x = "classLab", y = "log2CD16a", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("CD16a")), title = "CD16a_TPM", palette = "npg") +
      stat_compare_means(label =  "p.signif", comparisons = my_comparisons) + mytheme

ggboxplot(combineClass, x = "FcrnLab", y = "log2CD16a", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("CD16a")),title = "CD16a_TPM", palette = "npg") +
      stat_compare_means(label =  "p.signif", aes(group = IgG1CPKLab)) + mytheme

ggboxplot(combineClass, x = "classLab", y = "log2C1Q", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("C1Q")),title = "C1Q_TPM", palette = "npg") +
      stat_compare_means(label =  "p.signif", comparisons = my_comparisons) + mytheme

ggboxplot(combineClass, x = "FcrnLab", y = "log2C1Q", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("C1Q")),title = "C1Q_TPM", palette = "npg") +
      stat_compare_means(label =  "p.signif", aes(group = IgG1CPKLab)) + mytheme

dev.off()

##===============NK cell infiltration in FcRN and IgG1+IgG3 fraction============
combineClass$barcode = substr( rownames(combineClass), 1, 15)
FcRN_NKcell = merge(combineClass, NKcell, by.x = "barcode", by.y = 0)
colnames(FcRN_NKcell)[12:ncol(FcRN_NKcell)] = gsub(" |-|/","_", colnames(FcRN_NKcell)[12:ncol(FcRN_NKcell)])

pdf("Fig_NKcell_in_FcRN_and_IgG1fraction.allCancer.pdf", width = 3, height = 10)

FcRN_NKcell$log2NK_MCP = log2(FcRN_NKcell$NK_cell_MCP_COUNTER)
ggboxplot(FcRN_NKcell, x = "classLab", y = "log2NK_MCP", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("NK_cell")), title = "NKcell_MCP", palette = "npg") +
stat_compare_means(comparisons = my_comparisons,label =  "p.signif") + mytheme
ggboxplot(FcRN_NKcell, x = "FcrnLab", y = "log2NK_MCP", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("NK_cell")), title = "NKcell_MCP", palette = "npg") +
stat_compare_means(aes(group = IgG1CPKLab),label =  "p.signif") + mytheme

FcRN_NKcell$log2Macrophage_MCP = log2(FcRN_NKcell$Macrophage_Monocyte_MCP_COUNTER)
ggboxplot(FcRN_NKcell, x = "classLab", y = "log2Macrophage_MCP", ylim = c(1,14), width = 0.5, color = "IgG1CPKLab",ylab = expression(log[2]("Macrophage")), title = "Macrophage/Monocyte_MCP_COUNTER", palette = "npg") +
stat_compare_means(comparisons = my_comparisons,label =  "p.signif") + mytheme
ggboxplot(FcRN_NKcell, x = "FcrnLab", y = "log2Macrophage_MCP", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("Macrophage")), title = "Macrophage/Monocyte_MCP_COUNTER", palette = "npg") +
stat_compare_means(aes(group = IgG1CPKLab),label =  "p.signif") + mytheme

FcRN_NKcell$log2M2_xcell = log2(FcRN_NKcell$Macrophage_M2_XCELL)
ggboxplot(FcRN_NKcell, x = "classLab", y = "log2M2_xcell", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("Macrophage")), title = "M2cell_XCELL") +
stat_compare_means(comparisons = my_comparisons,label =  "p.signif") + mytheme

# Add by Li
#XCELL
FcRN_NKcell$log2NK_xcell = log2(FcRN_NKcell$NK_cell_XCELL)
ggboxplot(FcRN_NKcell, x = "classLab", y = "log2NK_xcell", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("NK_cell")), title = "XCELL") +
stat_compare_means(comparisons = my_comparisons,label =  "p.signif") + mytheme

FcRN_NKcell$log2M_xcell = log2(FcRN_NKcell$Macrophage_XCELL)
ggboxplot(FcRN_NKcell, x = "classLab", y = "log2M_xcell", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("Macrophage")), title = "XCELL") +
stat_compare_means(comparisons = my_comparisons,label =  "p.signif") + mytheme

# EPIC
FcRN_NKcell$tmp = log2(FcRN_NKcell$Macrophage_EPIC)
ggboxplot(FcRN_NKcell, x = "classLab", y = "tmp", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("Macrophage")), title = "EPIC") +
stat_compare_means(comparisons = my_comparisons,label =  "p.signif") + mytheme

FcRN_NKcell$tmp = log2(FcRN_NKcell$NK_cell_EPIC)
ggboxplot(FcRN_NKcell, x = "classLab", y = "tmp", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("NK_cell")), title = "EPIC") +
stat_compare_means(comparisons = my_comparisons,label =  "p.signif") + mytheme

# CIBERTSORT-ABS
FcRN_NKcell$tmp = log2(FcRN_NKcell$Macrophage_M0_CIBERSORT_ABS 
											 + FcRN_NKcell$Macrophage_M1_CIBERSORT_ABS
											 + FcRN_NKcell$Macrophage_M2_CIBERSORT_ABS)
ggboxplot(FcRN_NKcell, x = "classLab", y = "tmp", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("Macrophage")), title = "CIBERSORT-ABS") +
stat_compare_means(comparisons = my_comparisons,label =  "p.signif") + mytheme

FcRN_NKcell$tmp = log2(FcRN_NKcell$NK_cell_activated_CIBERSORT_ABS)
ggboxplot(FcRN_NKcell, x = "classLab", y = "tmp", width = 0.5, color = "IgG1CPKLab", ylab = expression(log[2]("activated_NK_cell")), title = "CIBERSORT-ABS") +
stat_compare_means(comparisons = my_comparisons,label =  "p.signif") + mytheme
dev.off()
