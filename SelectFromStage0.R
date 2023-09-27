#read in the seed amounts
Sinc<- read.csv("~/Documents/Wheat/2023/HarvestData/SI23circleHarvestFile.csv")
colnames(Sinc)[6]<- 'germplasmName'
Sinc<- Sinc[,c("germplasmName", "Weight..lb.")]

#read in notes
SI<- read.csv("~/Documents/GitHub/Wheat-Selection-Decisions-2023/SI_Urb_23-notes.csv")

#combine both
Sinc<- merge(Sinc, SI, by='germplasmName')
S0lines<- Sinc[c(grep('IL2022F4S', Sinc$germplasmName), grep('IL2022F3S', Sinc$germplasmName)),]
write.csv(S0lines, file='Stage0Selectionfile_2023.csv')
