setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2023")
library(asreml)
library(reshape)
library(doBy)

#read in means from stage 1 of analysis from previous years
previous<- read.csv('~/Documents/GitHub/Wheat-Selection-Decisions-2022/allILTrainingSetPhenotypes2022.csv', 
                    row.names=1)

#read in new trial results
data<- read.csv('YT blues 2023.csv', as.is=TRUE, row.names=1)
q<- matrix(unlist(strsplit(data$study, split="_")), nrow=3)
data$loc<- q[2,]
data$year<- q[3,]
data$site<- paste(data$loc, data$year, sep="_")
data$wt<- 1/(data$std.error^2)
data<- data[,colnames(previous)]

#combine both sources
all<- rbind(data, previous)

#standardize trait names
all[which(all$trait=="Grain.moisture.content.....CO_321.0001198"),'trait']<- 'Grain.moisture.content'
all[which(all$trait=="Heading.time...Julian.date..JD..CO_321.0001233"),'trait']<- 'Heading.time...Julian.date..JD..'
all[which(all$trait=="Lodging.incidence...0.9.percentage.scale.CO_321.0001282"),'trait']<- 'Lodging.incidence...0.9'
all[which(all$trait=="Maturity.time...spike.estimation...Julian.date..JD..CO_321.0501101"),'trait']<- 'Maturity.time'

#Multi-trait model with the relationship matrix
mod_test<- asreml(fixed=value.x~trait, 
            random=~trait:study+us(trait):germplasmName+trait:study:germplasmName,
            weights=value.y, data=allsub4, 
            na.action = na.method(y='include', x='include'), 
            family=asr_gaussian(dispersion = 1))
mod_test<- mkConv(mod_test)
blups_Agronomic_July14<- predict(mod_test, classify='germplasmName:trait')$pvals
write.csv(blups_Agronomic_July14, file='blups_Agronomic_July14.csv')


