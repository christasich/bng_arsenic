
library(pacman)
pacman::p_load(ggplot2,dplyr,rpart)


library(tidyverse)
library(boot)
library(leaps)
library(glmnet)

### LOAD DATA ###

setwd('C:/Users/chris/Documents/Projects/bng_arsenic/')

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
  select(As_ppb,depth=depth_ft,lon,lat,geocode) %>%
  mutate(As_10 = as.factor(ifelse(As_ppb>=10,1,0))) %>%
  mutate(As_50 = as.factor(ifelse(As_ppb>=50,1,0))) %>%
  mutate(as = ifelse(As_ppb==0,0.001,As_ppb),log_as=log(as))

d<-df
# load data
d <- read.csv("bamwsp.csv") %>%
  select(geocode, depth=depth_ft, as=arsenic_ppb) %>%
  mutate(as = ifelse(as==0, 0.001,as), # replace 0 with 0.001
         log_as = log(as)) # log of As 


# num of folds
k <- 10

# create index
folds <- rep_len(1:k, nrow(d))

# shuffle index using sample
folds <- sample(folds, nrow(d))

# number of models
nmodels <- 8

model.error <- matrix(data=NA, nrow=nrow(d), ncol=nmodels)
for (i in 1:k){
  # select rows in fold
  fold <- which(folds==i)
  
  # subset training data
  train <- d[-fold,]
  
  # model 1 = null model
  null <- mean(train$log_as)
  
  # model 2 = linear model
  linear <- lm(log_as ~ depth + lat + lon, data=train)
  
  # model 2.1 linear model
  linear2 <- lm(as ~ depth + lat + lon, data = train)
  
  # model 3 = regression tree
  rtree <- rpart(log_as ~ depth + lat + lon, data=train)
  
  # model 4 trees
  r2tree <- tree(log_as~ depth + lat + lon, data=train)
  
  # model 5 = random forest
  boost <- gbm(log_as ~ depth + lat + lon,data=train,distribution='gaussian',n.trees=1000,shrinkage=0.01,interaction.depth=4)
  
  # model 6 random forest
  rf <- randomForest(as~depth+lat+lon,data=train,ntree=501,nodesize=100)
  
  # test set
  test <- d[fold,]
  
  # store observed and predictions for folds
  model.error[fold,1] = exp(test$log_as)
  model.error[fold,2] = exp(null)
  model.error[fold,3] = exp(predict(linear,test))
  model.error[fold,4] = predict(linear2,test)
  model.error[fold,5] = exp(predict(rtree,test))
  model.error[fold,6] = exp(predict(r2tree,test))
  model.error[fold,7] = exp(predict(boost,test))
  model.error[fold,8] = predict(rf,test)
}

# store errors in a dataframe
model.results <- data.frame(obs=model.error[,1],
                            null=model.error[,2],
                            linear=model.error[,3],
                            linear2=model.error[,4],
                            rtree=model.error[,5],
                            r2tree=model.error[,6],
                            boost=model.error[,7],
                            rf=model.error[,8])

# rmse function
rmse <- function(obs,est){sqrt(mean((obs-est)^2,na.rm=T))}

# calculate rmse for models
rmse.all <- apply(model.results[,2:nmodels],2, rmse, obs=model.results[,1])

