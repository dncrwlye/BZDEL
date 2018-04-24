#......Alex all you should need to do is change the working directory
#setwd("/Users/buckcrowley/Desktop/BDEL/BZDEL/meta_analysis/")

setwd("C:/Users/r83c996/Documents/BZDEL/meta_analysis")

ncores = 4
tot.reps=500
reps.per.worker=round(tot.reps/ncores)
sampling_effort = TRUE
source('R/filo and hnv positive factorization/filo_factorization.R', echo=TRUE, print.eval=TRUE)

ncores = 4
tot.reps=500
reps.per.worker=round(tot.reps/ncores)
sampling_effort = TRUE
source('R/filo and hnv positive factorization/hnv.factorization.R')

ncores = 4
tot.reps=500
reps.per.worker=round(tot.reps/ncores)
sampling_effort = FALSE 
source('R/filo and hnv positive factorization/filo_factorization no sampling effort.R')

ncores = 4
tot.reps=500
reps.per.worker=round(tot.reps/ncores)
sampling_effort = FALSE 
source('R/filo and hnv positive factorization/hnv.factorization no sampling effort.R')



