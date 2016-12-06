### LOAD LIBRARIES ###

library(tidyverse)


### LOAD DATA ###

setwd('C:/Projects/Vanderbilt/bng_arsenic/')

as.df = tbl_df(read.csv('data/bamwsp.csv',header=T))

xy.df = tbl_df(read.csv('data/Mouza_25km.csv'))
xy.df = rename(xy.df,geocode = GEO2,lon = X_COORD,lat = Y_COORD)

df = left_join(as.df,xy.df,by='geocode')

df = df %>%
  filter(!is.na(lon)) %>%
  filter(!is.na(lat)) %>%
  select(arsenic_ppb,depth_ft,lon,lat) %>%
  mutate(as_10 = as.factor(ifelse(arsenic_ppb>=10,1,0))) %>%
  mutate(as_50 = as.factor(ifelse(arsenic_ppb>=50,1,0)))

summary(df['lon'])

### DIAGNOSTIC PLOTS OF DATA ###

# Pairs plot with 10% of data
df.sample = sample_frac(df,0.1,replace=F)
pairs(df.sample)

# Removed obs with 0 ppb as concentration to better see rest of data
df.sample1 = filter(df.sample,arsenic_ppb>0)
plot(df.sample1$arsenic_ppb~df.sample1$depth_ft)

# Histogram
df.sample2 = filter(df,arsenic_ppb>0)
hist(df.sample2$arsenic_ppb,breaks=60,xlab='As Concentration (ppb)',ylab='counts',col=2,main='Histogram of Non-Zero As Concentrations')

### PREDICT ARSENIC LEVEL OF NEW WELL ###

## LINEAR REGRESSION ##
set.seed(1)

df.train = sample_frac(df,0.7,replace=F)
df.test = setdiff(df,df.train)

lm.fit = lm(arsenic_ppb~.,data=df.train)
summary(lm.fit)

prob = predict(lm.fit,df.test,type='response')
pred.10 = ifelse(prob>=10,1,0)
pred.50 = ifelse(prob>=50,1,0)

cm.10 = table(df.test$as_10,pred.10)
cm.10

tp.10 = cm.10[4]/(cm.10[3]+cm.10[4])
tp.10
tn.10 = cm.10[1]/(cm.10[1]+cm.10[2])
tn.10

cm.50 = table(df.test$as_50,pred.50)
cm.50

tp.50 = cm.50[4]/(cm.50[3]+cm.50[4])
tp.50
tn.50 = cm.50[1]/(cm.50[1]+cm.50[2])
tn.50
