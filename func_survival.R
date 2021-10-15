## function for survival analysis
## install necessary libraries
p = c('survival','ggplot2')
for(el in p){
  if (!is.element(el, installed.packages()[,1]))install.packages(el, dep=TRUE)
  suppressWarnings(suppressMessages(invisible(require(el, character.only=TRUE))))
}

# define function
## uniform digit of sample id
getID<-function(sId,digits=3){
  tmp=sapply(sId,function(x){
    return(
      paste0(
        strsplit(x,split = '\\-')[[1]][1:digits],collapse='-')
    )})
  return(tmp)
}
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
getClass <- function(metaDat,first,cutoff=0.25){
  # define a cutoff
  low = quantile(metaDat[,first], cutoff , na.rm=T)
  high = quantile(metaDat[,first], 1-cutoff , na.rm=T)
  ### group sample into two groups
  metaDat$class = NA
  metaDat$class[ metaDat[,first] <= low ] = 1 # low fraction
  metaDat$class[ metaDat[,first] > high ] = 2 # high fraction
  metaDat = metaDat[ !is.na(metaDat$class),]
  return(metaDat)
}


## Collect data from bo and generate a datafrme to store required information for downstream analysis
generateData<-function(igRatio,clinDataPath,cancerTypePath,tumorPurityPath){
  ##########Input##########
  # load data
  tumorPurity=get(load(tumorPurityPath))
  clinData=get(load(clinDataPath))
  IgG1 = na.omit(igRatio)
  load(cancerTypePath)
  #cancerType = cancer.subtype
  cancerType = cancerType

  #########################
  #### uniform sample name
  ageBefore=clinData$age
  names(ageBefore)=clinData$patient

  survData=clin[,c(3,2)]
  rownames(survData)=clinData$patient

  names(IgG1)=getID(names(IgG1))
  names(cancerType)=getID(names(cancerType))
  names(tumorPurity)=getID(names(tumorPurity))

  ## Get overlapped samples
  tmp.ssAll=Reduce(intersect,list(rownames(survData[which(survData[,2]<=5000),]),
                                  names(IgG1),
                                  names(ageBefore),
                                  names(cancerType),
                                  names(tumorPurity)))
  # combine data
  result=data.frame(
    OS=survData[tmp.ssAll,2],
    Event=survData[tmp.ssAll,1],
    IgFraction=IgG1[tmp.ssAll],
    age=ageBefore[tmp.ssAll],
    cancer=as.character(cancerType[tmp.ssAll]),
    purity=tumorPurity[tmp.ssAll],
    stringsAsFactors=F)

  return(result)
}


figure6d<-function(metaDat,title, CUTOFF = 0.25){
  metaDat = getClass(metaDat,first='IgFraction',cutoff= CUTOFF)

  tmp.fit=survfit(Surv(metaDat$OS,metaDat$Event)~class,data=metaDat)
  plot(tmp.fit,col=c('darkblue','cornflowerblue'),main=title,lty=c(1,1),lwd=c(3,3),cex.lab=1.5,cex.axis=1.5,xlab='Days',ylab='Survival',mark.time=T,cex.main=1.5)
  legend('topright',legend=c('1:Low fraction','2:High fraction'),cex=0.8,lwd=4,col= c('darkblue','cornflowerblue'),bty='n')
  ## coxph regression, control age
  #1 vs 2
  case12=metaDat[metaDat$class %in% c(1,2),]
  cx12=coxph(Surv(case12$OS,case12$Event)~IgFraction+age+purity,data=case12)
  ## add coefficient to plot
  legend('bottomleft',legend=c(paste0('2 vs 1: HR=',signif(summary(cx12)$coef["IgFraction",'exp(coef)'],3),
                                      " P = ",
                                      signif(summary(cx12)$coef["IgFraction",'Pr(>|z|)'],3))),
         cex=0.8,text.col='red',bty='n')
}

getHarzardcoef<-function(metaDat,cutoff= 0.25, metric = "IgFraction"){
  ## for specific extration of coxph coef
  hazardAndpValue=c()

  getCoxphCoef<-function(sur_dat){
    surv=Surv(sur_dat$OS,sur_dat$Event)
    errflag=F
    sur_dat$class=as.character(sur_dat$class)

    if(length(unique(sur_dat$class))>1){
      if( unique(sur_dat$cancer) == "LAML"){
        cox.fit = tryCatch(coxph(surv~age+class, data=sur_dat),
                           error = function(e) errflag <<- T,
                           warning = function(w) errflag <<- F)
      }else{
        cox.fit = tryCatch(coxph(surv~age+class+purity, data=sur_dat),
                           error = function(e) errflag <<- T)
                           #warning = function(w) errflag <<- F) # supress warning for TGCT PDL1 and CTLA4
      }
      if(!errflag){
        regCoef = summary(cox.fit)$coef
        return(regCoef[grep('class',rownames(regCoef)),c('exp(coef)','Pr(>|z|)')])
      }
    }
    return(c(NA,NA))
  }

  for(cancer in unique(metaDat$cancer)){

    singleCancer = getClass(metaDat[metaDat$cancer==cancer,],first = metric,cutoff)

    ########
    if(any(table(singleCancer$class)<25))next
    #if(any(table(singleCancer$Event)<10))next
    print( cancer )
    case12 = singleCancer[singleCancer$class %in% c(1,2),]
    case12Coef = getCoxphCoef(case12)
    coxRegCoef = data.frame(Cancer=cancer,
                        Case='low vs high',
                        HR=case12Coef[1],
                        Pvalue=case12Coef[2],stringsAsFactors=F)
    if(!all(is.na(coxRegCoef$HR))){
      hazardAndpValue=rbind(hazardAndpValue,coxRegCoef)
    }
  }
  #hazardAndpValue$Pvalue = p.adjust(hazardAndpValue$Pvalue,method = 'BH')
  hazardAndpValue = hazardAndpValue[ order(hazardAndpValue$HR, decreasing=T),]
  return(hazardAndpValue)
}

getHarzardcoefWithTIB<-function(metaDat,cutoff= 0.25){
  ## for specific extration of coxph coef
  hazardAndpValue=c()

  getCoxphCoef<-function(sur_dat){
    surv=Surv(sur_dat$OS,sur_dat$Event)
    errflag=F
    sur_dat$class=as.character(sur_dat$class)

    if(length(unique(sur_dat$class))>1){
      cox.fit = tryCatch(coxph(surv~age+class+purity+Bcell, data=sur_dat),
                         error = function(e) errflag <<- T,
                         warning = function(w) errflag <<- T)
      if(!errflag){
        regCoef = summary(cox.fit)$coef
        return(regCoef[grep('class',rownames(regCoef)),c('exp(coef)','Pr(>|z|)')])
      }
    }
    return(c(NA,NA))
  }

  for(cancer in unique(metaDat$cancer)){

    singleCancer = getClass(metaDat[metaDat$cancer==cancer,],first='IgFraction',cutoff)

    ########
    if(any(table(singleCancer$class)<50))next
    #if(any(table(singleCancer$Event)<10))next
    case12 = singleCancer[singleCancer$class %in% c(1,2),]
    case12Coef = getCoxphCoef(case12)
    coxRegCoef = data.frame(Cancer=cancer,
                        Case='low fraction vs high fraction',
                        HR=case12Coef[1],
                        Pvalue=case12Coef[2],stringsAsFactors=F)
    if(!all(is.na(coxRegCoef$HR))){
      hazardAndpValue=rbind(hazardAndpValue,coxRegCoef)
    }
  }
  #hazardAndpValue$Pvalue = p.adjust(hazardAndpValue$Pvalue,method = 'BH')
  hazardAndpValue = hazardAndpValue[ order(hazardAndpValue$HR, decreasing=T),]
  return(hazardAndpValue)
}

getCxRes <- function(metaDat, metric = "IgFraction") {
  metaDat$IgFraction = metaDat[, metric]
  cxRes <- Reduce(rbind, lapply( unique(metaDat$cancer), function( x ){
    subdata = subset(metaDat, cancer == x)
    #print(nrow(subdata))
    if( nrow(subdata)<50 ){return(NA)}
    #subdata$IgFraction = subdata$IgFraction/max(subdata$IgFraction)
    if( x == "LAML" ){
      tmp.cx = coxph( Surv( OS, Event) ~ IgFraction+age, data = subdata)
    }else{
      subdata = subdata[ !is.na(subdata$purity), ]
      if( nrow(subdata)<50 ){ return(NA)}
      tmp.cx = coxph( Surv( OS, Event) ~ IgFraction+age+purity, data = subdata)
    }
    return( data.frame( Cancer = x ,HR = summary(tmp.cx)$coef["IgFraction","exp(coef)"],
                        Pvalue = summary(tmp.cx)$coef["IgFraction","Pr(>|z|)"]))
  }))
  return( na.omit(cxRes) )
}

getCxResRmTIB <- function(metaDat,rmMetric = "Bcell") {
  cxRes <- Reduce(rbind, lapply( unique(metaDat$cancer), function( x ){
    subdata = subset(metaDat, cancer == x)
    #print(nrow(subdata))
    if( nrow(subdata) < 50){return(NA)}
    #subdata$IgFraction = subdata$IgFraction/max(subdata$IgFraction)
    tmp.cx = coxph( Surv( OS, Event) ~ IgFraction+Bcell+age+purity, data = subdata)
    return( data.frame( Cancer = x ,HR = summary(tmp.cx)$coef["IgFraction","exp(coef)"],
                        Pvalue = summary(tmp.cx)$coef["IgFraction","Pr(>|z|)"]))
  }))
  return( na.omit(cxRes) )
}
