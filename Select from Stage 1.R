#read in the seed amounts
Sinc<- read.csv("~/Documents/Wheat/2023/HarvestData/SI23circleHarvestFile.csv")
colnames(Sinc)[6]<- 'germplasmName'

#read in notes
YT<- read.csv("~/Documents/GitHub/Wheat-Selection-Decisions-2023/YT_Urb_23-notes.csv")
SI<- read.csv("~/Documents/GitHub/Wheat-Selection-Decisions-2023/SI_Urb_23-notes.csv")

#read in gebvs
a<- read.csv('V2gebvs_Yield_July24.2023.csv')
a$trait<- 'Yield'
b<- read.csv('gebvs_Testweight_July24.2023.csv')
b$trait<-'TestWeight'
d<- read.csv('gebvs_Heading_July24.2023.csv')
d$trait<- 'Heading'
e<- read.csv('gebvs_Height_July24.2023.csv')
e$trait<- 'Height'
f<- read.csv('gebvs_Mat_July24.2023.csv')
f$trait<- 'Maturity'
g<- read.csv('gebvs_Lodg_July24.2023.csv')
g$trait<- 'Lodging'
h<- read.csv('GEBV_onlyDON.csv')[,c(1:5)]
h$trait<- 'DON'
all<- rbind(a, b, d, e, f, g, h)[,-1]

#subset lines of interest
cks<- c("07-19334","16-8048","Kaskaskia","Pio25R74")
lines<- all[grep('IL2022-', all$germplasmName),'germplasmName']
sub<- all[which(all$germplasmName %in% c(cks, lines)),]
sub$predicted.value<- round(sub$predicted.value, 2)

#make summary 
library(reshape)
smry<- cast(sub, germplasmName~trait, value='predicted.value')
smry2<- merge(smry, Sinc[,c(6,8)], by='germplasmName', all.x=TRUE)
smry3<- merge(smry2, SI[,c(1,3)], by='germplasmName', all.x=TRUE)
smry4<- merge(smry3, YT[,c(1,3)], by='germplasmName', all.x=TRUE)

write.csv(smry4, file='Stage1Selectionfile_2023.csv')

