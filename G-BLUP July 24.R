setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2023")
library(asreml)
library(reshape)
library(doBy)

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

#get data from 2020 to 2023
allsub<- droplevels.data.frame(all[which(all$year>19),])

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

save(allsub3, file='Blues2023.RData')

###Switch to the server from here
load('Blues2023.RData')

#trials to exclude because not int the TPE
ixel<- which(allsub3$study %in% c('Big6_Vin_23', 'YT_War_22', 'YT_War_23',  'Big6_Mas_23',  'YT_Msn_22'))
allsub3<- allsub3[-ixel,]

#####################################
#####    PREPARE MARKER DATA    #####
#####################################
library(asreml)
library(reshape)
library(gaston)
library(rrBLUP)
library(Matrix)

#Get the NEW marker data to include in analysis
setwd("/home/sharedFolder/2023_IL_GBS_Data")
geno_file23  <- "IL_all_regions_filt.vcf.gz"
geno23 <- read.vcf(geno_file23)

#Get the OLD data to include in analysis
setwd("/home/jrut")
geno_file  <- "IL_2022_all_regions_samp_filt_fullnames_dedup_imp.vcf.gz"
geno <- read.vcf(geno_file)
geno@ped$id<- sub("^.*:", "", geno@ped$id)

#correct the line names
yrseries<- gsub("20", "", as.character(2000:2019))
for(i in 1:length(yrseries)){
  geno@ped$id<- gsub(paste("IL", yrseries[i], sep=""), yrseries[i], geno@ped$id)
}
geno@ped$id<- gsub('IL20', "2020", geno@ped$id)
geno@ped$id<- gsub('IL21', "IL2021", geno@ped$id)
geno@ped$id<- gsub('16LCSDH', "IL16LCSDH", geno@ped$id)
geno@ped$id<-gsub("PIO-25R74", "Pio25R74", geno@ped$id)
geno@ped$id<-gsub("KASKASKIA", "Kaskaskia", geno@ped$id)

#combine both SNP sets
geno23sub<- select.snps(geno23, geno23@snps$id %in% geno@snps$id)
genosub<- select.snps(geno, geno@snps$id %in% geno23@snps$id)
genoAll<-rbind(genosub, geno23sub)

## subset the phenotypic data to include all those with genotypic data
inter_lines <- intersect(unique(allsub3$germplasmName), genoAll@ped$id)
allsub4<- droplevels.data.frame(allsub3[which(allsub3$germplasmName %in% inter_lines),])


## subset the genotypic data to include all those with phenotypic data + new lines
#training and validation lines
new_lines<- genoAll@ped$id[grep('IL2022-', genoAll@ped$id)]
genoAll2<- select.inds(genoAll, id %in% unique(c(inter_lines, new_lines)))

## subset the genotypic data to include only polymorphic snps
genoAll3<- select.snps(genoAll2, maf > 0.01)

#make relationship matrix
genoAll3<- as.matrix(genoAll3)-1
K<- A.mat(genoAll3)

#make K positive semidefinite 
K2<- nearPD(K)
K2<- K2$mat

#remove extra objects
rm(K)
rm(genoAll2)
rm(genoAll3)
rm(genoAll)
rm(genosub)
rm(geno23sub)

#####################################
#####        FIT MODELS         #####
#####################################

#save(allsub4, file='allsub4.RData')

#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}


#Multi-trait model without the relationship matrix
#mod_test<- asreml(fixed=value.x~trait, 
#            random=~trait:study+diag(trait):germplasmName,
#            weights=value.y, data=allsub4, 
#            na.action = na.method(y='include', x='include'), 
#            family=asr_gaussian(dispersion = 1))
#mod_test<- mkConv(mod_test)
#blups_Agronomic<- predict(mod_test, classify='germplasmName:trait')$pvals
#write.csv(blups_Agronomic, file='blups_Agronomic_July21.2023.csv')

#Multi-trait model with the relationship matrix
asreml.options(maxit= 150)
mod_gebv<- asreml(fixed=value.x~trait+trait:study, 
                  random=~diag(trait):vm(germplasmName, K2),
                  weights=value.y, data=allsub4, 
                  na.action = na.method(y='include', x='include'), 
                  family=asr_gaussian(dispersion = 1), workspace='64gb')
mod_gebv<- mkConv(mod_gebv)
gebvs_Agronomic<- predict(mod_gebv, classify='germplasmName:trait', pworkspace='8gb')$pvals
write.csv(gebvs_Agronomic, file='gebvs_Agronomic_July29.2023.csv')

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
                 random=~study+vm(germplasmName, K2),
                 weights=value.y, data=yld, 
                 na.action = na.method(y='include', x='include'), 
                 family=asr_gaussian(dispersion = 1), workspace='32gb')
#mod_yld2<- update(mod_yld, random.= ~.+germplasmName:study)
ldGEBV<- predict(mod_yld, classify='germplasmName', pworkspace='8gb')$pvals
write.csv(ldGEBV, file='V2gebvs_Yield_July24.2023.csv')

#Test Weight
mod_tw<- asreml(fixed=value.x~1, 
                random=~study+vm(germplasmName, K2),
                weights=value.y, data=tw, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
#mod_tw2<- update(mod_tw, random.= ~.+germplasmName:study)
twGEBV<- predict(mod_tw, classify='germplasmName', pworkspace='8gb')$pvals
write.csv(twGEBV, file='gebvs_Testweight_July24.2023.csv')

#Heading
mod_hd<- asreml(fixed=value.x~1, 
                random=~study+vm(germplasmName, K2),
                weights=value.y, data=hd, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
#mod_hd2<- update(mod_hd, random.= ~.+germplasmName:study)
hdGEBV<- predict(mod_hd, classify='germplasmName', pworkspace='8gb')$pvals
write.csv(hdGEBV, file='gebvs_Heading_July24.2023.csv')

#Height
#ht<- ht[-which(ht$value.x>500),]
mod_ht<- asreml(fixed=value.x~1, 
                random=~study+vm(germplasmName, K2),
                weights=value.y, data=ht, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
#mod_ht2<- update(mod_ht, random.= ~.+germplasmName:study)
htGEBV<- predict(mod_ht, classify='germplasmName', pworkspace='8gb')$pvals
write.csv(htGEBV, file='gebvs_Height_July24.2023.csv')

#Mat
mod_mat<- asreml(fixed=value.x~1, 
                random=~vm(germplasmName, K2),
                weights=value.y, data=mat, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
matGEBV<- predict(mod_mat, classify='germplasmName', pworkspace='8gb')$pvals
write.csv(matGEBV, file='gebvs_Mat_July24.2023.csv')


#Lodg
mod_lodg<- asreml(fixed=value.x~1, 
                random=~vm(germplasmName, K2),
                weights=value.y, data=lodg, 
                na.action = na.method(y='include', x='include'), 
                family=asr_gaussian(dispersion = 1), workspace='8gb')
lodgGEBV<- predict(mod_lodg, classify='germplasmName', pworkspace='8gb')$pvals
write.csv(lodgGEBV, file='gebvs_Lodg_July24.2023.csv')


