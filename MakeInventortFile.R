
#read in seed lots
SL<- read.csv('~/Documents/Wheat/2023/SI23_seedlotUpload.csv')

#stage 1 lines
stage1<- read.csv('~/Documents/GitHub/Wheat-Selection-Decisions-2023/Stage1Selectionfile_2023-selected.csv')
stage1selected<- stage1[which(stage1$decision==1),c('germplasmName')]

#stage 0 lines
stage0<- read.csv('~/Documents/GitHub/Wheat-Selection-Decisions-2023/Stage0Selectionfile_2023-selected.csv')
stage0selected<- stage0[which(stage0$Decision==1),c('germplasmName')]
allSelected<- c(stage1selected, stage0selected)

#add decision column
SL$Decision='DROP'
SL[which(SL$accession_name %in% allSelected),'Decision']<- 'KEEP'

head(SL)
SL<- SL[,c(1,2,5,9,10)]
pnumb<- as.numeric(gsub("SI_Urb_23-", "", SL$source))
SL<- SL[order(pnumb),]
SL$ent2324<- ""
SL[which(SL$Decision=='KEEP'),'ent2324']<- 1:length(which(SL$Decision=='KEEP'))
row.names(SL)<- 1:nrow(SL)

write.csv(SL, file='SI2023_SortFile.csv')

#read file
SIlist<- read.csv('~/Documents/Wheat/2023/SI2023_SortFile.csv')
SIlist<- SIlist[order(SIlist$ent2324),]
SIlist<- SIlist[which(SIlist$Decision=='KEEP'),]
group<- matrix(unlist(strsplit(SIlist$accession_name, split="-")), nrow=2)[1,]
SIlist$group<- group
SIlist<- SIlist[order(SIlist$group),]
RowInTray<- rep(c(1:20), 40)[1:782]
Tray<- sort(rep(1:40, 20))[1:782]
SIlist$RowInTray<- RowInTray
SIlist$Tray<- paste("PR", Tray, sep="")
write.csv(SIlist, file='~/Documents/Wheat/2023/PureRows_SIlist.csv')
