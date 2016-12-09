### LOAD LIBRARIES ###

library(tidyverse)
library(boot)
library(leaps)
library(glmnet)

### LOAD DATA ###

setwd('C:/Projects/Vanderbilt/bng_arsenic/')

as.df = tbl_df(read.csv('data/bamwsp_full.csv',header=T))
as.df = rename(as.df,As_ppb = arsenic_ppb)

xy.df = tbl_df(read.csv('data/Mouza_25km.csv'))
xy.df = rename(xy.df,geocode = GEO2,lon = X_COORD,lat = Y_COORD)

df = left_join(as.df,xy.df,by='geocode')
df$geocode = as.factor(as.df$geocode)

df = df %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(lat)) %>%
  filter(depth_ft>0) %>%
  filter(depth_ft<3000) %>%
  filter(As_ppb<1000) %>%
  select(As_ppb,depth_ft,lon,lat,geocode) %>%
  mutate(As_10 = as.factor(ifelse(As_ppb>=10,1,0))) %>%
  mutate(As_50 = as.factor(ifelse(As_ppb>=50,1,0)))

summary(df)

### VISUALIZE DATA ###

# Pairs plot with 10% of data
df.sample = sample_frac(df,0.1,replace=F)
pairs(~As_ppb+depth_ft+lon+lat,data=df.sample,upper.panel=NULL,pch='.')

# Histogram of >0 As Concentrations
df.gt0As = filter(df,As_ppb>0)
hist(df.gt0As$As_ppb,breaks=60,xlab='As Concentration (ppb)',ylab='counts',col=2,main='Histogram of Non-Zero As Concentrations')

## Plot As Concentrations vs Other Variables

# Plot As concentrations vs Depth
plot(df$As_ppb,-df$depth_ft,axes=F,ann=F,col='darkgrey')
axis(3)
axis(2)
box()
mtext('Arsenic Concentration (ppb)',side=3,line=3,cex = par("cex.lab"))
mtext('Depth (ft)',side=2,line=3,cex = par("cex.lab"))

# Plot As concentrations vs Lon
plot(df$lon,df$As_ppb,xlab = 'Longitude',ylab='Arsenic Concentration (ppb)',col='darkgrey')

# Plot As concentrations vs Lat
plot(df$lat,df$As_ppb,xlab = 'Latitude',ylab='Arsenic Concentration (ppb)',col='darkgrey')

### PREDICT ARSENIC LEVEL OF NEW WELL ###

## Determine variables using Best Subset Selection

regfit.full = regsubsets(As_ppb~depth_ft+depth_ft*lon+lon+poly(lat,3),df)
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

# Use all variables but errors are terrible

## LINEAR REGRESSION ##
set.seed(1)

lm.fit = glm(As_ppb~depth_ft+depth_ft*lon+lon+poly(lat,3),data=df)
cv.error = cv.glm(df,lm.fit,K=10)$delta[1]

fit = lm(As_ppb~depth_ft+depth_ft*lat+depth_ft*lon+lat*lon+lon+lat,data=df)

# Linear regression has a terrible fit

# Ridge Regression

x=model.matrix(As_ppb~depth_ft+lon+lat,df)[,-1]
y=df$As_ppb

grid = 10^seq(10,-2,length=100)
ridge.mod = glmnet(x,y,alpha=0,lambda=grid)
dim(coef(ridge.mod))
ridge.mod$lambda[50]


set.seed(1)
train = sample(1:nrow(x),nrow(x)/2)
test=(-train)
y.test=y[test]

ridge.mod=glmnet(x[train,],y[train],alpha=0,lambda=grid,thresh=1e-12)
ridge.pred=predict(ridge.mod,s=4,newx=x[test,])
mean((ridge.pred-y.test)^2)

set.seed(1)
cv.out=cv.glmnet(x[train,],y[train],alpha=0)
plot(cv.out)
bestlam=cv.out$lambda.min
bestlam

ridge.pred=predict(ridge.mod,s=bestlam,newx=x[test,])
mean((ridge.pred-y.test)^2)

out = glmnet(x,y,alpha=0)
predict(out,type='coefficients',s=bestlam)[1:4,]


# Lasso

lasso.mod = glmnet(x[train,],y[train],alpha=1,lambda=grid)
plot(lasso.mod)

set.seed(1)
cv.out=cv.glmnet(x[train,],y[train],alpha=1)
plot(cv.out)
bestlam=cv.out$lambda.min
lasso.pred=predict(lasso.mod,s=bestlam,newx=x[test,])
mean((lasso.pred-y.test)^2)

out = glmnet(x,y,alpha=1,lambda=grid)
lasso.coef=predict(out,type='coefficients',s=bestlam)[1:4,]
lasso.coef

# Classification

# Logistic Regression

log.fit = glm(As_ppb~depth_ft+lat+lon,data=df)
summary(log.fit)

cv.glm(df,log.fit,K=10)$delta[1]

#boost

library(gbm)

boost.as = gbm(As_ppb~depth_ft+lat+lon,distribution='gaussian',n.trees=100,shrinkage=0.01,interaction.depth=4,data=df)
