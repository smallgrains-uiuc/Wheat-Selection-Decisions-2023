setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2023")
library(asreml)
library(reshape)
library(gaston)
library(rrBLUP)
library(Matrix)

#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

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
genoAll2<- select.snps(genoAll, maf > 0.01)
genoAll3<- as.matrix(genoAll2)-1


#prepare phenotypic data
pheno<- read.csv("IL_Scb_data_for_analysis_2023.csv")
a<- read.csv('V2gebvs_Yield_July24.2023.csv')
pheno2<- pheno[which(pheno$germplasmName %in% c(a$germplasmName, c('Pio 25R74', '12-26004'))), ]
colnames(pheno2)[40]<- 'DON'
colnames(pheno2)[41]<- 'FDK'
pheno2$studyName<- as.factor(pheno2$studyName)
pheno2$germplasmName<- as.factor(pheno2$germplasmName)
pheno2$blockNumber<- as.factor(pheno2$blockNumber)

#make relationship matrix
genoScb<- genoAll3[c(a$germplasmName),]
write.csv(genoScb, file='genoScab23.csv')
#genoScb<- read.csv('genoScab23.csv')
K<- A.mat(genoScb)
K2<- nearPD(K)
K2<- K2$mat

#prepare phenotypic data
#pheno2<- pheno[which(pheno$germplasmName %in% genoScb[,1]),]
pheno2$studyName<- as.factor(pheno2$studyName)
pheno2$germplasmName<- as.factor(pheno2$germplasmName)
pheno2$blockNumber<- as.factor(pheno2$blockNumber)
#colnames(pheno2)[40]<- 'DON'
#colnames(pheno2)[41]<- 'FDK'
#row.names(K2)<- genoScb[,1]
#colnames(K2)<- genoScb[,1]
#########################################
### single-trait G-BLUP for scab resistance
#########################################
phenoDON<- pheno2[which(!is.na(pheno2$DON)),]
phenoDON[which(phenoDON$DON==-9),'DON']<- NA

amod<- asreml(fixed=DON~1,
              random= ~studyName+ 
              blockNumber:studyName+ 
              vm(germplasmName, K2),
              residual= ~dsum(~units | studyName),
              data=phenoDON,  
              na.action = na.method(y='include', x='include'), workspace="32gb")
GEBV_onlyDON<- predict(amod, classify='germplasmName', pworkspace="8gb")$pvals
pev<- GEBV_onlyDON[,'std.error']^2
Vg<- summary(amod)$varcomp['vm(germplasmName, K2)','component'] * mean(diag(K2))
GEBV_onlyDON$rel<- 1-(pev/Vg) 
write.csv(GEBV_onlyDON, file='GEBV_onlyDON.csv')

#########################################
### single-trait BLUP for scab resistance
#########################################

amod<- asreml(fixed=DON~1,
              random= ~studyName+ 
              blockNumber:studyName+ germplasmName,
              data=phenoDON,  
              na.action = na.method(y='include', x='include'), workspace="4gb")
BLUP_onlyDON<- predict(amod, classify='germplasmName', pworkspace="4gb")$pvals
pev<- BLUP_onlyDON[,'std.error']^2
Vg<- summary(amod)$varcomp['germplasmName','component'] 
BLUP_onlyDON$rel<- 1-(pev/Vg) 
write.csv(BLUP_onlyDON, file='BLUP_onlyDON.csv')



#check if predicted DON values correlate with this year's FDK data
pheno23<- pheno2[which(pheno2$studyYear==2023),]
pheno23_2<- merge(pheno23, GEBV_onlyDON, by='germplasmName')
q<- pheno23_2[intersect(which(pheno23_2$replicate=='1'), grep("19-", pheno23_2$germplasmName)),]
cor(q[,c('FDK', 'predicted.value','Heading.time...Julian.date..JD..CO_321.0001233')], use='pairwise.complete.obs')
q<- pheno23_2[intersect(which(pheno23_2$replicate=='2'), grep("19-", pheno23_2$germplasmName)),]
cor(q[,c('FDK', 'predicted.value','Heading.time...Julian.date..JD..CO_321.0001233')], use='pairwise.complete.obs')

#########################################
### multi-trait BLUP for scab resistance
#########################################
pheno2[which(pheno2$DON==-9),'DON']<- NA
pheno2[which(pheno2$FDK==-9),'FDK']<- NA

#Multi-trait model with FDK and DON
asreml.options(ai.sing=TRUE, aom=FALSE)
amod2<- asreml(fixed=cbind(DON, FDK)~trait,
               random= ~trait:studyName+ trait:studyName:blockNumber+us(trait):vm(germplasmName, K2),
               residual = ~id(units):us(trait),workspace=32e6, 
               data=pheno2,  na.action = na.method(y='include', x='include'))
summary(amod2)
GEBVscab<- predict(amod2, classify='germplasmName:trait', pworkspace="4gb")$pvals
write.csv(GEBVscab, file='GEBVscabJul17_2023.csv')

#get genetic correlations from model results
traits_sub<- c('DON', 'FDK')
gencor<- matrix(nrow=length(traits_sub), ncol=length(traits_sub))
for(i in 1:length(traits_sub)){
  for(j in 1:length(traits_sub)){
    num1<-summary(amod2)$varcomp[paste('trait:vm(germplasmName, K2)!trait_', traits_sub[i], ":", traits_sub[j], sep=""),'component']
    num2<-summary(amod2)$varcomp[paste('trait:vm(germplasmName, K2)!trait_', traits_sub[j], ":", traits_sub[i], sep=""),'component']
    num<- unique(na.omit(c(num1,num2)))
    denom_a<- summary(amod2)$varcomp[paste('trait:vm(germplasmName, K2)!trait_', traits_sub[i], ":", traits_sub[i], sep=""),'component']
    denom_b<- summary(amod2)$varcomp[paste('trait:vm(germplasmName, K2)!trait_', traits_sub[j], ":", traits_sub[j], sep=""),'component']
    cr<- num/c(sqrt(denom_a)*sqrt(denom_b))
    gencor[i,j]<- cr
  }
}
colnames(gencor)<- traits_sub
rownames(gencor)<- traits_sub

#matrix of genetic correlations
gencor

#calculate narrow sense heritability for each trait
ixG<- match(pheno2$germplasmName, row.names(K2))
meanD<- mean(diag(K2[ixG, ixG]))
Vg<- summary(amod2)$varcomp[paste('trait:vm(germplasmName, K2)!trait_', traits_sub, ":", traits_sub, sep=""),'component']
Ve<- summary(amod2)$varcomp[paste('units:trait!trait_', traits_sub, ":", traits_sub, sep=""),'component']
h2<- (meanD*Vg)/((meanD*Vg)+Ve)
names(h2)<- traits_sub
#vector of heritabilities
h2


##Check GS model accuracies by masking advanced lines

ix23<- which(pheno2$studyName == 'IL_Scb_23')
linestoMask<- unique(pheno2[ix23,][grep('19-', pheno2[ix23, 'germplasmName']),'germplasmName'])
phenoSim<- pheno2
phenoSim[which(phenoSim$germplasmName %in% linestoMask),'DON']<- NA

#Multi-trait model with FDK and DON
asreml.options(ai.sing=TRUE, aom=FALSE)
amodSim<- asreml(fixed=cbind(DON, FDK)~trait,
               random= ~trait:studyName+ trait:studyName:blockNumber+us(trait):vm(germplasmName, K2),
               residual = ~id(units):us(trait),workspace=32e6, 
               data=phenoSim,  na.action = na.method(y='include', x='include'))
summary(amodSim)
GEBVscabSim<- predict(amodSim, classify='germplasmName:trait', pworkspace="4gb")$pvals
GEBVscabSim<- GEBVscabSim[which(GEBVscabSim$germplasmName %in% linestoMask),]
GEBVscabSim<- GEBVscabSim[which(GEBVscabSim$trait=='DON'),]
trueDON<- GEBV_onlyDON[match(GEBVscabSim$germplasmName, GEBV_onlyDON$germplasmName),'predicted.value']
cor(GEBVscabSim$predicted.value, trueDON)

#Multi-trait model with FDK and DON, Exclude FDK data all years except this year
phenoSim2<- phenoSim
phenoSim2[which(phenoSim$germplasmName %in% linestoMask & phenoSim$blockNumber==1),'FDK']<- NA
phenoSim2[which(phenoSim$germplasmName %in% linestoMask & phenoSim$studyYear < 2023),'FDK']<- NA
asreml.options(ai.sing=TRUE, aom=FALSE)
amodSim<- asreml(fixed=cbind(DON, FDK)~trait,
                 random= ~trait:studyName+ trait:studyName:blockNumber+us(trait):vm(germplasmName, K2),
                 residual = ~id(units):us(trait),workspace=32e6, 
                 data=phenoSim2,  na.action = na.method(y='include', x='include'))
summary(amodSim)
GEBVscabSim<- predict(amodSim, classify='germplasmName:trait', pworkspace="4gb")$pvals
GEBVscabSim<- GEBVscabSim[which(GEBVscabSim$germplasmName %in% linestoMask),]
GEBVscabSim<- GEBVscabSim[which(GEBVscabSim$trait=='DON'),]
trueDON<- GEBV_onlyDON[match(GEBVscabSim$germplasmName, GEBV_onlyDON$germplasmName),'predicted.value']
cor(GEBVscabSim$predicted.value, trueDON)

