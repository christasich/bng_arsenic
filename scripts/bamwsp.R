### LOAD LIBRARIES ###
library(pacman)
pacman::p_load(ggplot2,dplyr,rpart)
library(e1071)
library(leaps)

### LOAD DATA ###

setwd('C:/Projects/Vanderbilt/bng_arsenic/')

as.df = tbl_df(read.csv('data/bamwsp.csv',header=T))
as.df = rename(as.df,As_ppb = arsenic_ppb)

xy.df = tbl_df(read.csv('data/Mouza_25km.csv'))
xy.df = rename(xy.df,geocode = GEO2,lon = X_COORD,lat = Y_COORD)

df = left_join(as.df,xy.df,by='geocode')
df$geocode = as.factor(as.df$geocode)

d = df %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(lat)) %>%
  filter(depth_ft>0) %>%
  filter(depth_ft<3000) %>%
  filter(As_ppb<1000) %>%
  select(As_ppb,depth=depth_ft,lon,lat,geocode) %>%
  mutate(As_10 = as.factor(ifelse(As_ppb>=10,1,0))) %>%
  mutate(As_50 = as.factor(ifelse(As_ppb>=50,1,0))) %>%
  mutate(as = ifelse(As_ppb==0,0.001,As_ppb),log_as=log(as))

### VISUALIZE DATA ###

# Pairs plot with 10% of data
d.sample = sample_frac(d,0.1,replace=F)
pairs(~as+log_as+lon+lat,data=d.sample,upper.panel=NULL,pch='.')

# Histogram of >0 As Concentrations
d.gt0As = filter(d,as>0)
hist(d.gt0As$as,breaks=60,xlab='As Concentration (ppb)',ylab='counts',col=2,main='Histogram of Non-Zero As Concentrations')

## Plot As Concentrations vs Other Variables

# Plot As concentrations vs Depth
plot(d$as,-d$depth,axes=F,ann=F,col='darkgrey')
axis(3)
axis(2)
box()
mtext('Arsenic Concentration (ppb)',side=3,line=3,cex = par("cex.lab"))
mtext('Depth (ft)',side=2,line=3,cex = par("cex.lab"))

# Plot As concentrations vs Lon
plot(d$lon,d$as,xlab = 'Longitude',ylab='Arsenic Concentration (ppb)',col='darkgrey')

# Plot As concentrations vs Lat
plot(d$as,d$lat,ylab = 'Latitude',xlab='Arsenic Concentration (ppb)',col='darkgrey')

## Determine variables using Best Subset Selection

regfit.full = regsubsets(as~depth+lon+lat,d)
reg.summary = summary(regfit.full)

# Plot model statistics

par(mfrow=c(2,2))
plot(reg.summary$rss ,xlab="Number of Variables",ylab="RSS",type="l")
best.rss = which.min(reg.summary$rss)
points(best.rss,reg.summary$rss[best.rss],col="red",cex=2,pch=20)

plot(reg.summary$adjr2 ,xlab =" Number of Variables ",ylab=" Adjusted RSq",type="l")
best.r2 = which.max(reg.summary$adjr2)
points(best.r2,reg.summary$adjr2[best.r2],col="red",cex=2,pch=20)

plot(reg.summary$cp ,xlab ="Number of Variables",ylab="Cp",type='l')
best.cp=which.min(reg.summary$cp)
points (best.cp,reg.summary$cp[best.cp],col="red",cex=2,pch=20)

plot(reg.summary$bic,xlab ="Number of Variables",ylab="BIC",type='l')
best.bic = which.min(reg.summary$bic)
points(best.bic,reg.summary$bic[best.bic],col="red",cex=2,pch=20)

par(mfrow=c(2,2))
plot(regfit.full,scale ="r2")
plot(regfit.full,scale ="adjr2")
plot(regfit.full,scale ="Cp")
plot(regfit.full,scale ="bic")

### MODELS ###

# num of folds
k <- 10

# create index
folds <- rep_len(1:k, nrow(d))

# shuffle index using sample
folds <- sample(folds, nrow(d))

# number of models
nmodels <- 11

model.error <- matrix(data=NA, nrow=nrow(d), ncol=nmodels)
for (i in 1:k){
  # select rows in fold
  fold <- which(folds==i)
  
  # subset training data
  train <- d[-fold,]
  
  # model 1 = null model of As
  null1 <- mean(train$as)
  
  # model 2 = null model of log(As)
  null2 <- mean(train$log_as)
  
  # model 3 = linear model of As
  linear1 <- lm(as ~ depth + lat + lon, data=train)
  
  # model 4 = linear model of log(As)
  linear2 <- lm(log_as ~ depth + lat + lon, data=train)
  
  # model 5 = regression tree of As
  rtree1 <- rpart(as ~ depth + lat + lon, data=train)
  
  # model 6 = regression tree of log(As)
  rtree2 <- rpart(log_as ~ depth + lat + lon, data=train)
  
  # model 7 = polynomial model
  poly1 <- lm(as ~ poly(depth,3) + lat + lon, data=train)
  
  # model 8 = polynomial model
  poly2 <- lm(as ~ poly(depth,4) + lat + lon, data=train)

  # model 9 = polynomial model
  poly3 <- lm(as ~ depth + poly(lat,2) + lon, data=train)
  
  # model 10 = SVM
  svm <- svm(as ~ depth + lat + lon, data=train)
  
  # test set
  test <- d[fold,]
  
  # store observed and predictions for folds
  model.error[fold,1] = exp(test$log_as)
  model.error[fold,2] = null1
  model.error[fold,3] = exp(null2)
  model.error[fold,4] = predict(linear1,test)
  model.error[fold,5] = exp(predict(linear2,test))
  model.error[fold,6] = predict(rtree1,test)
  model.error[fold,7] = exp(predict(rtree2,test))
  model.error[fold,8] = predict(poly1,test)
  model.error[fold,9] = predict(poly2,test)
  model.error[fold,10] = predict(poly3,test)
  model.error[fold,11] = predict(svm,test)
}

# store errors in a dataframe
model.results <- data.frame(obs=model.error[,1],
                            null1=model.error[,2],
                            null2=model.error[,3],
                            linear1=model.error[,4],
                            linear2=model.error[,5],
                            rtree1=model.error[,6],
                            rtree2=model.error[,7],
                            poly1=model.error[,8],
                            poly2=model.error[,9],
                            poly3=model.error[,10],
                            svm=model.error[,11])

# rmse function
rmse <- function(obs,est){sqrt(mean((obs-est)^2,na.rm=T))}

# calculate rmse for models
rmse.all <- apply(model.results[,2:nmodels],2, rmse, obs=model.results[,1])

### PLOTS ###

plot(rmse.all,xaxt='n',xlab='Model',ylab='RMSE',ylim=range(0:max(rmse.all)))
axis(1,at=1:10,labels=names(rmse.all))
points(5,rmse.all[5],col="red",cex=2,pch=20)

# Best Model

plot(rtree1, uniform=TRUE, main="Regression Tree for Arsenic Concentration")
text(rtree1, all=TRUE, cex=.8)
