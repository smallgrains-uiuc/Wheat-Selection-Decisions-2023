setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2023")
library(asreml)
library(reshape)
library(doBy)

###########Analysis is to select among stage 3 and 4 lines (Advanced)

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


#subset the the Advanced (stage 3 and 4 lines) and checks
cks<- c("07-19334", "12-26004", "16-8048", "Kaskaskia", 'Pio 25R74', 'AgriMAXX503','Pio 25R76')
yt<- read.csv('YT blues 2023.csv', as.is=TRUE, row.names=1)
Adv<- yt[which(yt$study=='YT_Blv_23'),'germplasmName']
allsub4<- droplevels.data.frame(allsub3[which(allsub3$germplasmName %in% unique(c(cks, Adv))),])
allsub4$group<- 4
allsub4[grep("19-", allsub4$germplasmName),'group']<- 3
allsub4[which(allsub4$germplasmName %in% cks),'group']<- 0

#trials to include
ixel<- which(allsub4$study %in% c('Big6_Vin_23', 'YT_Msn_22'))
allsub4<- droplevels.data.frame(allsub4[-ixel,])

#Fit each trait separately
yld<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Grain.yield...bu.ac'),]))
tw<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Grain.test.weight...lbs.bu'),]))
hd<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Heading.time...Julian.date..JD..'),]))
ht<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Plant.height.inches'),]))
lodg<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Lodging.incidence...0.9'),]))
mat<- na.omit(droplevels.data.frame(allsub4[which(allsub4$trait== 'Maturity.time'),]))

#Yield
asreml.options(maxit= 150)
mod_yld<- asreml(fixed=value.x~1, 
                 random=~study+fa(study, 2):germplasmName,
                 weights=value.y, data=yld, 
                 na.action = na.method(y='include', x='include'), 
                 family=asr_gaussian(dispersion = 1), workspace='4gb')
mod_yld<- mkConv(mod_yld)
yldBLUP<- predict(mod_yld, classify='study:germplasmName', pworkspace='4gb')$pvals
yldBLUP$trait<- 'Yield'
yldBLUPavg<- predict(mod_yld, classify='germplasmName', average='study', pworkspace='4gb')$pvals
yldBLUPavg$trait<- 'Yield'
yldBLUPavg$study<- 'Average'

#Test Weight
mod_tw<- asreml(fixed=value.x~1, 
                random=~study+fa(study, 2):germplasmName,
                weights=value.y, data=tw, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
#mod_tw2<- update(mod_tw, random.= ~.+germplasmName:study)
twBLUP<- predict(mod_tw, classify='study:germplasmName', pworkspace='8gb')$pvals
twBLUP$trait<- 'TestWeight'
twBLUPavg<- predict(mod_tw, classify='germplasmName', average='study',pworkspace='4gb')$pvals
twBLUPavg$trait<- 'TestWeight'
twBLUPavg$study<- 'Average'

#Heading
mod_hd<- asreml(fixed=value.x~1, 
                random=~study+germplasmName,
                weights=value.y, data=hd, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
#mod_hd2<- update(mod_hd, random.= ~.+germplasmName:study)
hdBLUP<- predict(mod_hd, classify='germplasmName', pworkspace='8gb')$pvals
hdBLUP$trait<- 'Heading'
hdBLUP$study<- NA
hdBLUP<- hdBLUP[,colnames(twBLUP)]

#Height
mod_ht<- asreml(fixed=value.x~1, 
                random=~study+germplasmName,
                weights=value.y, data=ht, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
htBLUP<- predict(mod_ht, classify='germplasmName', pworkspace='8gb')$pvals
htBLUP$trait<- 'Height'
htBLUP$study= NA
htBLUP<- htBLUP[,colnames(twBLUP)]

#Mat
colnames(mat)[7]<- "predicted.value"
colnames(mat)[8]<- "std.error"
mat$status=NA
mat$study=NA
matBLUP<- mat[,colnames(htBLUP)]

#Lodg
colnames(lodg)[7]<- "predicted.value"
colnames(lodg)[8]<- "std.error"
lodg$status=NA
lodgBLUP<- lodg[,colnames(matBLUP)]

#DON
don<- read.csv("BLUP_onlyDON.csv")
don$trait<- 'DON'
don$study<- NA
donBLUP<- don[,colnames(matBLUP)]
donBLUP<- donBLUP[which(donBLUP$germplasmName %in% c(cks, Adv)),]

#add the SMBV data
sbmv<- read.csv("SBMVphenotype_download.csv")
colnames(sbmv)[40]<- 'SBMV'
sbmv$studyName<- as.factor(sbmv$studyName)
sbmv$germplasmName<- as.factor(sbmv$germplasmName)
mod_sbmv<- asreml(fixed=SBMV~1, 
                random=~studyName+germplasmName, data=sbmv, 
                na.action = na.method(y='include', x='include'),
                workspace='4gb')
sbmvBLUP<- predict(mod_sbmv, classify='germplasmName', pworkspace='8gb')$pvals
sbmvBLUP$trait<- 'SBMV'
sbmvBLUP$study<- NA


#combine all the BLUPS
all<- rbind(yldBLUPavg, twBLUPavg, htBLUP,hdBLUP, matBLUP, lodgBLUP, donBLUP, sbmvBLUP, yldBLUP, twBLUP)
all$predicted.value<- round(all$predicted.value, 2)

#make wide format
head(allsub4)
all_wide<- cast(all, germplasmName~trait+study, value='predicted.value')



#Add notes to the summary
#get notes
notes<- read.csv('notes_phenotype_download.csv')
notes<- notes[-grep('harvest day', notes$notes),]
notes$notes<- gsub(" (Operator: jrutkoski, Time: )", "", notes$notes, fixed=TRUE)
notes$notes<- trimws(notes$notes)
notes<- notes[-which(notes$notes==''),]
study<- gsub("_", "", notes$studyName)
notes$notes<- paste(study, paste("[", notes$notes, "]", sep=""), sep="")
ugid<- unique(notes$germplasmName)
gidnote<-c()
for(i in 1:length(ugid)){
  gidnote<- append(gidnote, paste(notes[which(ugid[i]== notes$germplasmName),'notes'], collapse="; "))
}
df<- data.frame(germplasmName=ugid, notes=gidnote)
all_wide<- merge(all_wide, df, by='germplasmName', all.x=TRUE, all.y=FALSE)
all_wide<- all_wide[,c(1,31,9,3,7,4,2,5,6,8, 32:52,10:30,53)]
write.csv(all_wide, file='SelectFromStage4_2023.csv')
