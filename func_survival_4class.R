p = c('survival','ggplot2')
for(el in p){
  if (!is.element(el, installed.packages()[,1]))install.packages(el, dep=TRUE)
  suppressWarnings(suppressMessages(invisible(require(el, character.only=TRUE))))
}

get2Class <- function(metaDat,first,second,cutoff=3){
  # define a cutoff
  shm.Cut = fivenum(metaDat[,first])

  Ig.Cut_low = fivenum(metaDat[metaDat[,first] <= shm.Cut[cutoff],second])
  Ig.Cut_high = fivenum(metaDat[metaDat[,first] > shm.Cut[6-cutoff],second])
  ### group sample into two groups
  metaDat$class = NA
  ## lOW SHM
  metaDat$class[metaDat[,first] <= shm.Cut[cutoff] & metaDat[,second] <= Ig.Cut_low[cutoff]] = 1
  metaDat$class[metaDat[,first] <= shm.Cut[cutoff]  & metaDat[,second] > Ig.Cut_low[6-cutoff]] = 2
  ## High SHM
  metaDat$class[metaDat[,first] > shm.Cut[6-cutoff] & metaDat[,second] <= Ig.Cut_high[cutoff]] = 3
  metaDat$class[metaDat[,first] > shm.Cut[6-cutoff] & metaDat[,second] > Ig.Cut_high[6-cutoff]] = 4
  metaDat = metaDat[!is.na(metaDat$class),]
  return(metaDat)
}

get2ClassBy37 <- function(metaDat,first,second,cutoff=0.3){
  # define a cutoff
  low = quantile(metaDat[,first], cutoff , na.rm=T)
  high = quantile(metaDat[,first], 1-cutoff , na.rm=T)

  Ig.low_low = quantile(metaDat[metaDat[,first] <= low,second], cutoff , na.rm=T)
  Ig.low_high = quantile(metaDat[metaDat[,first] <= low,second], 1-cutoff , na.rm=T)
  Ig.high_low = quantile(metaDat[metaDat[,first] > high,second], cutoff , na.rm=T)
  Ig.high_high = quantile(metaDat[metaDat[,first] > high,second], 1-cutoff , na.rm=T)

  ### group sample into two groups
  metaDat$class = NA
  ## Low SHM
  metaDat$class[metaDat[,first] <= low & metaDat[,second] <= Ig.low_low] = 1
  metaDat$class[metaDat[,first] <= low  & metaDat[,second] > Ig.low_high] = 2
  ## High SHM
  metaDat$class[metaDat[,first] > high & metaDat[,second] <= Ig.high_low] = 3
  metaDat$class[metaDat[,first] > high & metaDat[,second] > Ig.high_high] = 4
  metaDat = metaDat[!is.na(metaDat$class),]
  return(metaDat)
}

## Collect data from bo and generate a datafrme to store required information for downstream analysis
generate2Data<-function(shm,clinDataPath,igFraction,cancerTypePath,tumorPurityPath){
  ##########Input##########
  # load data
  tumorPurity=get(load(tumorPurityPath))
  clinData=get(load(clinDataPath))
  load(cancerTypePath)
  cancerType = cancer.subtype
  ### Uniform Name

  ageBefore=clinData$age
  names(ageBefore)=clinData$patient
  survData=clin[,c(3,2)]
  rownames(survData)=clinData$patient

  cancerTypeName=getID(names(cancerType))
  cancerType=as.character(cancerType)
  names(cancerType)=cancerTypeName

  shmName=getID(names(shm))
  shmData=as.numeric(shm)
  names(shmData)=shmName
  shmData=shmData[!is.na(shmData)]

  names(tumorPurity)=getID(names(tumorPurity))
  names(igFraction)=getID(names(igFraction))

  ## Get Overlap Sample
  tmp.ss=intersect(names(igFraction),names(shmData))
  tmp.ssAll=Reduce(intersect,list(rownames(survData[which(survData[,2]<=5000),]),
                                  tmp.ss,
                                  names(ageBefore),
                                  names(cancerType),
                                  names(tumorPurity)))
  # combine data
  result=data.frame(
    OS=survData[tmp.ssAll,2],
    Event=survData[tmp.ssAll,1],
    SHM=shmData[tmp.ssAll],
    Switch=igFraction[tmp.ssAll],
    age=ageBefore[tmp.ssAll],
    cancer=as.character(cancerType[tmp.ssAll]),
    purity=tumorPurity[tmp.ssAll],
    stringsAsFactors=F)

  return(result)
}

figure5c<-function(metaDat,cutoff=3,title){


  metaDat = get2Class(metaDat,first='SHM',second = 'Switch', cutoff = cutoff)

  tmp.fit=survfit(Surv(metaDat$OS,metaDat$Event)~class,data=metaDat)
  plot(tmp.fit,col=c('pink','lightblue','darkred','cornflowerblue'),lty=c(1,1,1,1),lwd=c(3,3,3,3),cex.lab=1.5,cex.axis=1.5,xlab='Days',ylab='Survival',mark.time=T,cex.main=1.5,main=title)
  legend('topright',legend=c('1:Low expr + Low fraction','2:Low expr + High fraction','3:High expr +Low fraction', '4:High expr +High fraction'),cex=0.8,lwd=4,col= c('pink','lightblue','darkred','cornflowerblue'),bty='n')
  ## coxph regression, control age
  #1 vs 2,lOW SHM
  case12=metaDat[metaDat$class %in% c(1,2),]
  cx12=coxph(Surv(case12$OS,case12$Event)~class+age+purity,data=case12)
  #3 vs 4,High SHM
  case34=metaDat[metaDat$class %in% c(3,4),]
  cx34=coxph(Surv(case34$OS,case34$Event)~class+age+purity,data=case34)
  ## add coefficient to plot
  legend('bottomleft',legend=c(paste0('2 vs 1: HR=',signif(summary(cx12)$coef["class",'exp(coef)'],3),
                                      ", p = ",
                                      signif(summary(cx12)$coef["class",'Pr(>|z|)'],3)),
                               paste0('4 vs 3: HR=',signif(summary(cx34)$coef["class",'exp(coef)'],3),
                                      ", p = ",
                                      signif(summary(cx34)$coef["class",'Pr(>|z|)'],3))),
         cex=0.8,text.col='red',bty='n')
}
