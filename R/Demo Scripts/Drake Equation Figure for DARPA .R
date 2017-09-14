library(tidyverse)


vector0 <- rep(c('barrier'), each=20)

vector1 <- c(rnorm(n= 10, mean= 0, sd = 10), rnorm(n = 10, mean=0, sd =2 ))


dataframe <- matrix(c(vector0, vector1), ncol=2)
                    
dataframe1 <- as.data.frame(dataframe)  %>%
  mutate(V2 = as.character(V2)) %>%
  mutate(V2= as.numeric(V2))

ggplot(dataframe1, aes(factor(V1), V2)) +
  geom_violin(aes(fill = V1))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        axis.ticks.x = element_blank(), axis.text = element_blank()) +
  coord_flip() +
  xlab('variance') +
  ylab('barrier to spillover') 

