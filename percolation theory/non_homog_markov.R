#.........................nipah markov chain simulator........................................................
library(tidyverse)
#.........................theoritical data....................................................................

#temperature 
day <- seq(1:365)
pi_x <- seq(from = 0, to= pi, by = .01)

temp_year <- sin(day/(365/2)*pi)+10 + rnorm(n=365,mean=0, sd=.3)
plot(temp_year, type='l')

shedding_year = temp_year * .4 + rnorm(n=365, mean=0, sd=.5)

plot(shedding_year)
plot(shedding_year, temp_year)
summary(lm(shedding_year~temp_year) )

persistence_year = exp(1/temp_year)  + rnorm(n=365, mean=0,sd=.01) 
plot(persistence_year, temp_year)
summary(lm(persistence_year~temp_year))

human_contact_year = dnorm(day, mean = 240, sd=40)  + rnorm(n=365, mean=0,sd=.001) 
plot(human_contact_year)

dataset <- cbind(day, temp_year, shedding_year, persistence_year, human_contact_year)
dataset <- as.data.frame(dataset) %>%
  mutate(shedding_cat = ifelse(shedding_year > 4,1,0)) %>%
  mutate(human_contact_year_cat = ifelse(human_contact_year > .005,1,0)) %>%
  mutate(persistence_year_cat = ifelse(persistence_year > 1.11,1,0))

persistence<-glm(persistence_year_cat ~ temp_year, data= dataset, family=binomial (link=logit)) %>%summary()
shedding <- glm(shedding_cat ~ temp_year, data= dataset, family=binomial (link=logit)) %>%summary()
human_contact <- glm(human_contact_year_cat ~ temp_year, data= dataset, family=binomial (link=logit)) %>%summary()


#...........................markov_chain............................................
starting_vector <- c(100,0,0,0)

temp_year <- c(temp_year, temp_year, temp_year, temp_year, temp_year)

storage_vector <- matrix(0,4,365*5)
for (i in 1:365*5)
{
  temp_day_i = temp_year[i]
  
  #..........................log odds ................................................
  
  shedding_p_i      = exp(shedding$coefficients[2,1] * temp_day_i)/ (1+ exp(shedding$coefficients[2,1]* temp_day_i))
  persistence_p_i   = exp(persistence$coefficients[2,1] * temp_day_i)/ (1+ exp(persistence$coefficients[2,1]* temp_day_i))
  human_contact_p_i = exp(human_contact$coefficients[2,1] * temp_day_i)/ (1+ exp(human_contact$coefficients[2,1]* temp_day_i))
  
  markov_matrix <-matrix(c(1-shedding_p_i, shedding_p_i, 0, 0, 
                           0, persistence_p_i, human_contact_p_i, 1- (human_contact_p_i + persistence_p_i),
                           0,0,1,0,
                           0,0,0,1), nrow=4, byrow = FALSE)
  
  starting_vector <- markov_matrix %*% (starting_vector)
  starting_vector[1,1] <- 100
  storage_vector[,i] <- starting_vector
}

storage_vector <- as.data.frame(t(storage_vector))

plot(storage_vector$V3)


lines(storage_vector$V2, type='l', col='red')
lines(storage_vector$V3, type='l', col='blue')
lines(storage_vector$V4, type='l', col='gren')


