p = c('ggplot2','scales')
for(el in p){
  if (!is.element(el, installed.packages()[,1]))install.packages(el, dep=TRUE)
  suppressWarnings(suppressMessages(invisible(require(el, character.only=TRUE))))
}

mytheme <- theme_light()+theme(axis.title = element_text(size=13,face = 'bold'),
                               axis.text = element_text(size=10,face = 'bold'),
                               axis.text.x = element_text(size=10,face = 'bold',angle=-90,vjust = 0.5),
                               legend.text = element_text(size=7,face = 'bold'),
                               panel.grid.major = element_blank(),
                               plot.title = element_text(size=15,face='bold',hjust = 0.5))

mytheme_noXlab <- theme(legend.title = element_text(size = 7),
                 legend.text = element_text(size=5),
                 legend.key.width=unit(1, "lines"),
                 plot.title = element_text(size=10, face="bold", hjust=0.5),
                 plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),
                 axis.text.x = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.y = element_text(face="bold"),
                 panel.background = element_rect(fill = "white", colour = "black"),
                 panel.grid.major = element_line(linetype = "solid",colour = "grey90"))

mytheme_grid = theme(legend.title = element_text(size = 7),
                 legend.text = element_text(size=5, face = "bold"),
                 legend.key.width=unit(1, "lines"),
                 plot.title = element_text(size=10, face="bold", hjust=0.5),
                 plot.margin = unit(c(1, 5, 0.5, 0.5), "lines"),
                 axis.text.x = element_text(angle = 90, size=9,face="bold",color="black"),
                 axis.title.x=element_blank(),
                 axis.text.y = element_text(face="bold"),
                 panel.background = element_rect(fill = "white", colour = "black"),
                 panel.grid.major = element_line(linetype = "solid",colour = "grey90"))


## Data format for <dat>:
#------------------------------------
# Cancer   Case         HR     Pvalue
# GBM 1 vs 2 2.93227849 0.55567483
# LUAD 1 vs 2 0.29340367 0.29789292
# COAD 1 vs 2 0.09920154 0.02622474
# COAD 3 vs 4 0.14332852 0.59383129
#------------------------------------
# x_axis (character): column NAME in <dat> you want to put on x axis
# y_axis (character): column NAME in <dat> you want to put on y axis
# rho (character): column NAME in <dat> your want to fill in each box
# mean_value (numeric): mean value of your heatmap color bar
# min_value (numeric): min value of your heatmap color bar
# max_value (numeric): max value of your heatmap color bar
# bar_anno (vector): values you want to show on your heatmap color bar
# box_size (numeric): size of box in your heatmap
# fdr_cutoff (numeric): cutoff of significant fdr

corPlot <- function(dat,x_axis,y_axis,rho,pvalue,fdr_cutoff = 0.1,
                    min_value=-1,mean_value=0,max_value=1,
                    bar_anno=c(-1,0,1),titleName='NULL',
                    box_size=9,barwidth = 0.8, barheight = 5, theme = mytheme){


  dat$sig = NA
  dat$sig[dat[,pvalue] > fdr_cutoff] = "FDRnotSig"
  dat$sig[dat[,pvalue] <= fdr_cutoff] = "FDRsig"
  dat$sig = factor(dat$sig,levels = c( "FDRnotSig", "FDRsig"))
  dat = dat[ !is.na(dat[,pvalue]) & !is.na(dat[,rho]),]

  ggplot(dat, aes_string(x_axis, y_axis )) +
    geom_point(aes_string(color = rho, size= 'sig' ), shape=15) +
    scale_size_manual(name='',values = c( "FDRnotSig" = 5, "FDRsig" = 9), label = c( paste0("FDR >= ",fdr_cutoff),paste0("FDR < ",fdr_cutoff)), drop = FALSE)+
    coord_equal()+
    labs(x='',y='',title=titleName)+
    scale_color_gradientn(colours = c("navy","white","firebrick3"),
                          values = rescale(c(min_value,mean_value,max_value)),
                          guide = "colorbar", limits=c(min_value,max_value),breaks=bar_anno,na.value="white")+
    scale_fill_gradientn(colours = c("navy","white","firebrick3"),
                         values = rescale(c(min_value,mean_value,max_value)),
                         guide = "colorbar", limits=c(min_value,max_value),breaks=bar_anno,na.value="white")+
    theme+guides(fill = guide_colorbar(barwidth = barwidth, barheight = barheight),
                  color = guide_colorbar(barwidth = barwidth, barheight = barheight))
}
