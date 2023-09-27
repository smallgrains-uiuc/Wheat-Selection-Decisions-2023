setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/Wheat-Selection-Decisions-2023")
library(asreml)
library(reshape)
library(doBy)

#read in new Big6 trial results
data<- read.csv('Big6 blues Aug29.2023.csv', as.is=TRUE, row.names=1)
q<- matrix(unlist(strsplit(data$study, split="_")), nrow=3)
data$loc<- q[2,]
data$year<- q[3,]
data$site<- paste(data$loc, data$year, sep="_")
data$wt<- 1/(data$std.error^2)

#make summary table
data$predicted.value2<- round(data$predicted.value, 2)
wide<- cast(data, germplasmName~trait+study, value='predicted.value2')
write.csv(wide, file='Big6 blues Aug29.2023-wide.csv')

######YIELD data 
dataGY<- droplevels.data.frame(data[which(data$trait=='Grain.yield...bu.ac'),])
#combined model without the relationship matrix
dataGY$site<- as.factor(dataGY$site)
dataGY$germplasmName<- as.factor(dataGY$germplasmName)
dataGY$wt<- c(1/(dataGY$std.error^2))
dataGY<- na.omit(dataGY)
mod<- asreml(fixed=predicted.value~1, 
             random=~site+fa(site, 2):germplasmName,
             weights=wt, 
             data=dataGY, 
             na.action = na.method(y='include', x='include'), 
             family=asr_gaussian(dispersion = 1))
mod<- mkConv(mod)
yldBLUP.Big6<- predict(mod, classify='site:germplasmName')$pvals
yldBLUP.Big6$predicted.value2<- round(yldBLUP.Big6$predicted.value,2)
wide.faBLUP<- cast(yldBLUP.Big6, germplasmName~site, value='predicted.value2')
write.csv(wide.faBLUP, file='Big6.Yield-FABLUPs.csv')

######Days to Heading data 
dataDH<- droplevels.data.frame(data[which(data$trait=='Heading.time...Julian.date'),])
#combined model without the relationship matrix
dataDH$site<- as.factor(dataDH$site)
dataDH$germplasmName<- as.factor(dataDH$germplasmName)
dataDH$wt<- c(1/(dataDH$std.error^2))
dataDH<- na.omit(dataDH)
mod<- asreml(fixed=predicted.value~1, 
             random=~site+fa(site, 2):germplasmName,
             weights=wt, 
             data=dataDH, 
             na.action = na.method(y='include', x='include'), 
             family=asr_gaussian(dispersion = 1))
mod<- mkConv(mod)
dhBLUP.Big6<- predict(mod, classify='site:germplasmName')$pvals
dhBLUP.Big6$predicted.value2<- round(dhBLUP.Big6$predicted.value,2)
wide.faBLUPdh<- cast(dhBLUP.Big6, germplasmName~site, value='predicted.value2')
wide.faBLUPdh<- wide.faBLUPdh[,c("germplasmName", 'Urb_23', 
                 'Prn_23', 'Woo_23', 'Mas_23', 'Laf_23', 'Vin_23')]
write.csv(wide.faBLUPdh, file='Big6.Heading-FABLUPs.csv')

######Heigth data 
dataHT<- droplevels.data.frame(data[which(data$trait=='Plant.height.inches'),])
#combined model without the relationship matrix
dataHT$site<- as.factor(dataHT$site)
dataHT$germplasmName<- as.factor(dataHT$germplasmName)
dataHT$wt<- c(1/(dataHT$std.error^2))
dataHT<- na.omit(dataHT)
mod<- asreml(fixed=predicted.value~1, 
             random=~site+fa(site, 2):germplasmName,
             weights=wt, 
             data=dataHT, 
             na.action = na.method(y='include', x='include'), 
             family=asr_gaussian(dispersion = 1))
mod<- mkConv(mod)
htBLUP.Big6<- predict(mod, classify='site:germplasmName')$pvals
htBLUP.Big6$predicted.value2<- round(htBLUP.Big6$predicted.value,2)
wide.faBLUPht<- cast(htBLUP.Big6, germplasmName~site, value='predicted.value2')
wide.faBLUPht<- wide.faBLUPht[,c("germplasmName", 'Urb_23', 
                                 'Prn_23', 'Woo_23', 'Mas_23')]
write.csv(wide.faBLUPht, file='Big6.Height-FABLUPs.csv')


######Test weigtwt data 
datatwt<- droplevels.data.frame(data[which(data$trait=='Grain.test.weight...lbs.bu'),])
#combined model without the relationship matrix
datatwt$site<- as.factor(datatwt$site)
datatwt$germplasmName<- as.factor(datatwt$germplasmName)
datatwt$wt<- c(1/(datatwt$std.error^2))
datatwt<- na.omit(datatwt)
mod<- asreml(fixed=predicted.value~1, 
             random=~site+fa(site, 2):germplasmName,
             weights=wt, 
             data=datatwt, 
             na.action = na.method(y='include', x='include'), 
             family=asr_gaussian(dispersion = 1))
mod<- mkConv(mod)
twtBLUP.Big6<- predict(mod, classify='site:germplasmName')$pvals
twtBLUP.Big6$predicted.value2<- round(twtBLUP.Big6$predicted.value,2)
wide.faBLUPtwt<- cast(twtBLUP.Big6, germplasmName~site, value='predicted.value2')
wide.faBLUPtwt<- wide.faBLUPtwt[,c("germplasmName", "Fre_23", 'Urb_23', "Fra_23",
                                   'Prn_23', 'Woo_23', 'Mas_23', 'Laf_23', 'Vin_23')]
write.csv(wide.faBLUPtwt, file='Big6.Testweight-FABLUPs.csv')


#Site correlations
#site correlation matrix
met.corr <-function(object,site,faN=2){ # ,faRS=1
  
  n<-nlevels(site)
  
  varcomp<-summary(object)$varcomp['component']
  vcn<-row.names(varcomp)
  aimn<-vcn[grep('fa\\(.*,.*\\)',vcn)]
  varcomp1<-varcomp[aimn,]
  vect1<-varcomp1[1:n]
  w.var<-diag(vect1)
  vect2<-varcomp1[(n+1):((1+faN)*n)]
  L.var<-matrix(vect2,nrow=n)
  wL.var<-L.var%*%t(L.var)+w.var
  df<-wL.var
  for(i in 1:(n-1)){
    for(j in 2:n){
      if(i<j){df[i,j]<-df[j,i]/(sqrt(df[i,i]*df[j,j]))
      j<-j+1}
    }
    i<-i+1
  }
  rownames(df)<-levels(site)
  colnames(df)<-levels(site)
  
  df.2<-df
  
  for(i in 1:(n-1)){
    for(j in 2:n){
      if(i<j){df[j,i]<-df[i,j]
      j<-j+1}
    }
    i<-i+1
  }
  diag(df)<-1
  return(df)
}
df<- met.corr(mod, unique(dataGY$site, 2))
heatmap(df, symm=TRUE)
df<- round(df, 2)
write.csv(df, 'siteCorrelations.Big6-2023.csv')

