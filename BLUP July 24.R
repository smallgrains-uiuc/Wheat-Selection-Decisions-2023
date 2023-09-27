setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2023")
library(asreml)
library(reshape)
library(doBy)

###########Analysis is to select among stage 3 and 4 lines


#########################################
#####    PREPARE PHENOTYPIC DATA    #####
#########################################

#read in means from stage 1 of analysis from previous years
previous<- read.csv('~/Documents/GitHub/Wheat-Selection-Decisions-2022/allILTrainingSetPhenotypes2022.csv', row.names=1)

#read in new YT trial results
data<- read.csv('YT blues 2023.csv', as.is=TRUE, row.names=1)
q<- matrix(unlist(strsplit(data$study, split="_")), nrow=3)
data$loc<- q[2,]
data$year<- q[3,]
data$site<- paste(data$loc, data$year, sep="_")
data$wt<- 1/(data$std.error^2)
data<- data[,colnames(previous)]

#combine both sources
all<- rbind(data, previous)

#read in new Big6 trial results
data<- read.csv('Big6 blues June24.2023.csv', as.is=TRUE, row.names=1)
q<- matrix(unlist(strsplit(data$study, split="_")), nrow=3)
data$loc<- q[2,]
data$year<- q[3,]
data$site<- paste(data$loc, data$year, sep="_")
data$wt<- 1/(data$std.error^2)
data<- data[,colnames(previous)]

#combine both sources
all<- rbind(all, data)

#standardize trait names
all[which(all$trait=="Grain.moisture.content.....CO_321.0001198"),'trait']<- 'Grain.moisture.content'
all[which(all$trait=="Heading.time...Julian.date..JD..CO_321.0001233"),'trait']<- 'Heading.time...Julian.date..JD..'
all[which(all$trait=="Lodging.incidence...0.9.percentage.scale.CO_321.0001282"),'trait']<- 'Lodging.incidence...0.9'
all[which(all$trait=="Maturity.time...spike.estimation...Julian.date..JD..CO_321.0501101"),'trait']<- 'Maturity.time'

#get data from 2021 to 2023
allsub<- droplevels.data.frame(all[which(all$year>20),])

#select traits to use
traits_sub<- c('Grain.yield...bu.ac', 'Grain.test.weight...lbs.bu',
               'Plant.height.inches','Heading.time...Julian.date..JD..',
               'Lodging.incidence...0.9','Maturity.time')
allsub2<- droplevels.data.frame(allsub[which(allsub$trait %in% traits_sub),])

#subset the yield trials
if(length(grep('Scb', unique(allsub2$study)))>0){
  trials_sub<- unique(allsub2$study)[-grep('Scb', unique(allsub2$study))]
  allsub2<- droplevels.data.frame(allsub2[which(allsub2$study %in% trials_sub),])
}

#cast and melt data to include missing data
q<- melt(cast(allsub2, site+loc+year+study+germplasmName~trait, value='predicted.value'), 
         variable_name='trait', 
         id.vars=c('site','loc','year','study','germplasmName'), na.rm=FALSE)
q2<- melt(cast(allsub2, site+loc+year+study+germplasmName~trait, value='wt'), 
          variable_name='trait', 
          id.vars=c('site','loc','year','study','germplasmName'), na.rm=FALSE)
allsub3<- merge(q, q2, by=c('site','loc','year','study','germplasmName', 'trait'), 
                all.x=TRUE, all.y=TRUE)

#convert variables to factors
ix<- which(colnames(allsub3) %in%c('site','loc','year','study','germplasmName', 'trait'))
for(i in ix){
  allsub3[,i ]<- as.factor(as.character(allsub3[, i]))
}

#change NA weights to 0
if(length(which(is.na(allsub3$value.y)))>0){
  allsub3[which(is.na(allsub3$value.y)),'value.y']<- 0
}


#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}


#subset the advanced lines only
Adv_lines<- unique(allsub3[which(allsub3$study=='YT_Blv_23'),'germplasmName'])
allsub4<- allsub3[which(allsub3$germplasmName %in% Adv_lines),]

#remove trials with no experimental lines
allsub4<- allsub4[-which(allsub4$study %in% c('Adv_Urb_21', 'Adv_Neo_21', 'Adv_Stj_21','YT_Msn_22','Big6_Vin_23')),]

#Fit each trait separately
yld<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Grain.yield...bu.ac'),]))
tw<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Grain.test.weight...lbs.bu'),]))
hd<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Heading.time...Julian.date..JD..'),]))
ht<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Plant.height.inches'),]))
lodg<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Lodging.incidence...0.9'),]))
mat<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Maturity.time'),]))


#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}


#Yield
asreml.options(maxit= 150)
mod_yld<- asreml(fixed=value.x~1, 
                 random=~study+fa(study, 2):germplasmName,
                 weights=value.y, data=yld, 
                 na.action = na.method(y='include', x='include'), 
                 family=asr_gaussian(dispersion = 1), workspace='4gb')
mod_yld<- mkConv(mod_yld)
ldGEBV<- predict(mod_yld, classify='study:germplasmName', pworkspace='4gb')$pvals
write.csv(ldGEBV, file='FAblups_Yield_July24.2023.csv')

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
df<- met.corr(mod_yld, unique(yld$study, 2))
heatmap(df)
df<- round(df, 2)
write.csv(df, 'siteCorrelations.ILadvanced-2023.csv')

#Test Weight
mod_tw<- asreml(fixed=value.x~1, 
                random=~study+fa(study, 2):germplasmName,
                weights=value.y, data=tw, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='4gb')
mod_tw<- mkConv(mod_tw)
twGEBV<- predict(mod_tw, classify='study:germplasmName', pworkspace='4gb')$pvals
write.csv(twGEBV, file='FAblups_Testweight_July24.2023.csv')

#Heading
mod_hd<- asreml(fixed=value.x~1, 
                random=~study+fa(study, 2):germplasmName,
                weights=value.y, data=hd, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='4gb')
mod_hd<- mkConv(mod_hd)
hdGEBV<- predict(mod_hd, classify='study:germplasmName', pworkspace='4gb')$pvals
write.csv(hdGEBV, file='FAblups_Heading_July24.2023.csv')

#Height
mod_ht<- asreml(fixed=value.x~1, 
                random=~study+fa(study, 2):germplasmName,
                weights=value.y, data=ht, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='4gb')
mod_ht<- mkConv(mod_ht)
htGEBV<- predict(mod_ht, classify='study:germplasmName', pworkspace='4gb')$pvals
write.csv(htGEBV, file='FAblups_Height_July24.2023.csv')

#Mat
mod_mat<- asreml(fixed=value.x~1, 
                 random=~germplasmName,
                 weights=value.y, data=mat, 
                 na.action = na.method(y='include', x='include'), 
                 family=asr_gaussian(dispersion = 1), workspace='4gb')
matGEBV<- predict(mod_mat, classify='germplasmName', pworkspace='4gb')$pvals
write.csv(matGEBV, file='blups_Mat_July24.2023.csv')


#Lodg
mod_lodg<- asreml(fixed=value.x~1, 
                  random=~germplasmName,
                  weights=value.y, data=lodg, 
                  na.action = na.method(y='include', x='include'), 
                  family=asr_gaussian(dispersion = 1), workspace='4gb')
lodgGEBV<- predict(mod_lodg, classify='germplasmName', pworkspace='4gb')$pvals
write.csv(lodgGEBV, file='blups_Lodg_July24.2023.csv')



