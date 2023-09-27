setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2023")

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

#make summary 
library(reshape)
smry<- cast(all, germplasmName~trait, value='predicted.value')

#Net merit calculation
wheat_price0<- mean(c(9.9128, 7.0402, 5.4621, 4.9414, 4.9757, 4.4014, 4.3945))
soybean_price<- mean(c(16.1773, 13.6890, 9.5344, 8.9298, 9.3456, 9.7820, 9.8753))
#wheat price fcn
wheatPrice<- function(don, twt, wheat_price0){
  donDiscount<- don*-0.2
  twtDiscount<- c(58-twt)*-.08
  wheat_price<- wheat_price0+donDiscount+twtDiscount
  return(wheat_price)
}
#net merit function
netMerit<- function(headings, yields, dons,  twts, wheat_price0, soybean_price){
  wheat_price1<- wheatPrice(dons, twts, wheat_price0)
  soy_yld_gain<- 0.5* (135-headings)
  soy_profit_gain<- soy_yld_gain*soybean_price
  wheat_profit<- yields*wheat_price1
  total_profit<- wheat_profit + soy_profit_gain
  return(total_profit)
}
merit<- netMerit(smry$Heading, smry$Yield, smry$DON, smry$TestWeight,wheat_price0, soybean_price)
smry$merit<- merit
cor(smry$merit, smry$TestWeight, use='pairwise.complete.obs')
cor(smry$merit, smry$Yield, use='pairwise.complete.obs')
cor(smry$merit, smry$Heading, use='pairwise.complete.obs')
cor(smry$merit, smry$DON, use='pairwise.complete.obs')

#subset lines that are still being tested
blue<- read.csv('YT blues 2023.csv', as.is=TRUE, row.names=1)
smry$InTest<- smry$germplasmName %in% blue$germplasmName

write.csv(smry, file='2023 parent selection file.csv')
