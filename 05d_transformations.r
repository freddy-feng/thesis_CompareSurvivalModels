rm(list = ls())
library(msigno)
set.cf()

library(moments)
load('data/05a_datasets.RData')

df_transf = df_long
p = nrow(long_vars)

for (i in 1:p) {
  varname = long_vars$variable[i]
  sk = skewness(df_long[ , varname], na.rm = T)
  long_vars$skew[i] = sk
  title = paste(i, ') ', varname, ', skewness = ', round(sk, 2), sep = '')
  hist(df_long[ , varname], 20, 
       main = title)
}

ids_skewpos = c(1:2, 4:5, 7:20, 22, 37:38)
# all these variables are non-negative :)

par(mfrow = c(2, 1)) # view these vertically!

for (i in ids_skewpos) {
  varname = long_vars$variable[i]
  df_transf[ , varname] = log(df_long[ , varname] + 1)
  hist(df_long[ , varname], 20, main = paste(i, ') ', varname))
  hist(df_transf[ , varname], 20, main = 'log(x+1)')
}

ids_skewneg = c(29:30, 31:32, 36)

for (i in 29:30) {
  varname = long_vars$variable[i]
  df_transf[ , varname] = df_long[ , varname]^3
  hist(df_long[ , varname], 20, main = paste(i, ') ', varname))
  hist(df_transf[ , varname], 20, main = 'x^3')
}

min(df_long[ , long_vars$variable[31]], na.rm = T)
min(df_long[ , long_vars$variable[32]], na.rm = T)
min(df_long[ , long_vars$variable[36]], na.rm = T)

for (i in 31:32) {
  varname = long_vars$variable[i]
  df_transf[ , varname] = (35+df_long[ , varname])^3
  hist(df_long[ , varname], 20, main = paste(i, ') ', varname))
  hist(df_transf[ , varname], 20, main = '(35+x)^3')
}

for (i in 36) {
  varname = long_vars$variable[i]
  df_transf[ , varname] = (1700+df_long[ , varname])^3
  hist(df_long[ , varname], 20, main = paste(i, ') ', varname))
  hist(df_transf[ , varname], 20, main = '(1700+x)^3')
}

# normalize everything
for (i in 1:p) {
  varname = long_vars$variable[i]
  df_transf[ , varname] = df_transf[ , varname] / sd(df_transf[ , varname], na.rm = T)
  long_vars$skew_new[i] = skewness(df_transf[ , varname], na.rm = T)
}

dfround(long_vars[ids_skewpos, c(1, 5:6)], 2)
dfround(long_vars[ids_skewneg, c(1, 5:6)], 2)

rm.but(c('df_transf', 'df_surv', 'long_vars'))
save.image('data/05d_datasets_transf.RData')



