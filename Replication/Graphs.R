
library(stringr)
library(caret)

#IMPORT DATA
expo_mac = read.csv("Leximin_Expomatch.csv")
gm_mac = read.csv("gm_mac.csv")
supr_expo_mac = read.csv("supr_Expo_mac.csv")


#CREATE DISTANCE COLUMN AND REMOVE FAKE VALUES
expo_mac = expo_mac[which(expo_mac$time.taken != 0),]
expo_mac$distance = expo_mac$Ground.truth.ATE - expo_mac$Estimated.ATE

gm_mac = gm_mac[which(gm_mac$time.taken != 0),]
gm_mac$distance = gm_mac$Ground.truth.ATE - gm_mac$Estimated.ATE

supr_expo_mac = supr_expo_mac[which(supr_expo_mac$time.taken != 0),]
supr_expo_mac$distance = supr_expo_mac$Ground.truth.ATE - supr_expo_mac$Estimated.ATE

#CLASSIFYING THE DATASETS INTO CATEGORIES
gm_mac[,1] = gsub("[[:digit:]]+", '', gm_mac$name)
gm_mac[,1] = gsub("-", '', gm_mac[,1])
gm_mac[,1] = gsub("\\.", '', gm_mac[,1])

expo_mac[,1] = gsub("[[:digit:]]+", '', expo_mac$name)
expo_mac[,1] = gsub("-", '', expo_mac[,1])
expo_mac[,1] = gsub("\\.", '', expo_mac[,1])

supr_expo_mac[,1] = gsub("[[:digit:]]+", '', supr_expo_mac$name)
supr_expo_mac[,1] = gsub("-", '', supr_expo_mac[,1])
supr_expo_mac[,1] = gsub("\\.", '', supr_expo_mac[,1])


unique(leximin_Expomatch[,1])

df = data.frame(matrix(ncol=16,nrow=0, dimnames=list(NULL, c("name", "Mean TE GenMatch", "Mean TE ExpoMatch",
                                                             "Mean TE Alt_ExpoMatch","lower Confint GenMatch","upper Confint GenMatch",
                                                             'lower Confint ExpoMatch','upper Confint ExpoMatch',
                                                             "lower Confint Alt_ExpoMatch","upper Confint Alt_ExpoMatch",
                                                             "lower Range of Pvals GenMatch","upper Range of Pvals GenMatch",
                                                             "lower Range of Pvals ExpoMatch","upper Range of Pvals ExpoMatch",
                                                             "lower Range of Pvals Alt_ExpoMatch","upper Range of Pvals Alt_ExpoMatch"
                                                             ))))
for (i in unique(expo_mac[,1])){
  exp_subset = expo_mac[which(expo_mac[,1] == i),]
  supr_subset = expo_mac[which(supr_expo_mac[,1] == i),]
  gm_subset = expo_mac[which(gm_mac[,1] == i),]
  
  print(i)
  print(quantile(gm_subset$distance,c(0.025,0.975)))
  
  df[nrow(df)+1,] = c(i,mean(gm_subset$distance),mean(exp_subset$distance),mean(supr_subset$distance),
                      quantile(gm_subset$distance,c(0.025,0.975)),quantile(exp_subset$distance,c(0.025,0.975)),quantile(supr_subset$distance,c(0.025,0.975)),
                      range(gm_subset$Best.p.value),range(exp_subset$Best.p.value),range(supr_subset$Best.p.value))
}

write.csv(df,"mac.csv")






layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
plot(gm_mac$Best.p.value,gm_mac$distance,
     main='p-value optimized value with distance to treatment effect',
     ylab='Discrepancy to treatment effect',
     xlab='Largest smallest p-value',ylim=c(-0.5,0.5),col="red")
points(supr_expo_mac$Best.p.value,supr_expo_mac$distance,
       main='p-value optimized value with distance to treatment effect',
       ylab='Discrepancy to treatment effect',
       xlab='Largest smallest p-value',ylim=c(-0.5,0.5),col="blue")
points(expo_mac$Best.p.value,expo_mac$distance,
     main='p-value optimized value with distance to treatment effect',
     ylab='Discrepancy to treatment effect',
     xlab='Largest smallest p-value',ylim=c(-0.5,0.5),col="green")
legend(0.67,0.6, legend=c("GenMatch", "Standard Expomatch","Alternative ExpoMatch"), 
       fill = c("red","blue","green"),cex=0.8)
plot(gm_mac$Supremum.value,gm_mac$distance,
     main='Supremum optimized value with distance to treatment effect',
     ylab='Discrepancy to treatment effect',
     xlab='Supremum value',ylim=c(-0.5,0.5),xlim=c(0,100),col="red")
points(supr_expo_mac$Supremum.value,supr_expo_mac$distance,
       main='p-value optimized value with distance to treatment effect',
       ylab='Discrepancy to treatment effect',
       xlab='Largest smallest p-value',ylim=c(-0.5,0.5),col="blue")
points(expo_mac$Supremum.value,expo_mac$distance,
      main='p-value optimized value with distance to treatment effect',
      ylab='Discrepancy to treatment effect',
      xlab='Largest smallest p-value',ylim=c(-0.5,0.5),col="green")
legend(72,0.6, legend=c("GenMatch", "Standard Expomatch","Alternative ExpoMatch"), 
       fill = c("red","blue","green"),cex=0.8)

expo_mac$name[which(expo_mac$Supremum.value > 100)]

supr_expo_mac$name[which(supr_expo_mac$Supremum.value > 100)]
gm_mac$name[which(gm_mac$Supremum.value > 100)]

ks.test(gm_mac$Best.p.value,supr_expo_mac$Best.p.value)
ks.test(expo_mac$Best.p.value,supr_expo_mac$Best.p.value)
t.test(gm_mac$Best.p.value,supr_expo_mac$Best.p.value)
t.test(expo_mac$Best.p.value,supr_expo_mac$Best.p.value)


ks.test(gm_mac$Supremum.value,supr_expo_mac$Supremum.value)
ks.test(expo_mac$Supremum.value,supr_expo_mac$Supremum.value)
t.test(gm_mac$Supremum.value,supr_expo_mac$Supremum.value)
t.test(expo_mac$Supremum.value,supr_expo_mac$Supremum.value)

ks.test(gm_mac$distance,supr_expo_mac$distance)
ks.test(expo_mac$distance,supr_expo_mac$distance)
t.test(gm_mac$distance,supr_expo_mac$distance)
t.test(expo_mac$distance,supr_expo_mac$distance)

plot(expo_mac$time.taken,supr_expo_mac$time.taken)

dim(expo_mac)
dim(supr_expo_mac)

hist(expo_mac$time.taken)
hist(supr_expo_mac$time.taken)
