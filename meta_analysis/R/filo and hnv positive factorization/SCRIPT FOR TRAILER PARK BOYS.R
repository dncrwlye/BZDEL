#......Alex all you should need to do is change the working directory
setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")
ncores = 7
tot.reps=500
reps.per.worker=round(tot.reps/ncores)

source('R/filo and hnv positive factorization/filo_factorization.R', echo=TRUE, print.eval=TRUE)

ncores = 7
tot.reps=500
reps.per.worker=round(tot.reps/ncores)

source('R/filo and hnv positive factorization/filo_factorization no sampling effort.R')

ncores = 7
tot.reps=500
reps.per.worker=round(tot.reps/ncores)

source('R/filo and hnv positive factorization/hnv_factorization.R')

ncores = 7
tot.reps=500
reps.per.worker=round(tot.reps/ncores)

source('R/filo and hnv positive factorization/hnv_factorization no sampling effort.R')



