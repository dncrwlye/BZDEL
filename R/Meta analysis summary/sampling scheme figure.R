#.................. sampling scheme................
library(ggplot2)
library(magrittr)
library(dplyr)
x<-seq(0,365,by=1)
amp <- rnorm(n=length(x), mean=.2)
amp <- runif(n=length(x))
y<-amp*sin(.04*x)+amp
y<-y*sd(y)
plot(y, type='l')

smooth_y <- vector()
start = 20
exp_value =1.1
for (i in seq(start,length(y)))
{
l <- i-(start-1)
vector_exp <- 1/(exp_value^(2:(start+1)))
smooth_y[i] <- sum(vector_exp * y[i:l]) * (1/sum(vector_exp))
}
plot(y, type='l')
lines(smooth_y, col='red', lwd=5)

smooth_yt <- as.data.frame((as.matrix(smooth_y))) %>%
  mutate(x = row_number())

rect_sampling_point_unique <- c(30, 30*3,  30*5,  30*7,  30*9, 30*11)
rect_sampling_point_unique <- data.frame(
  xmin = rect_sampling_point_unique,
  xmax = rect_sampling_point_unique + 8,
  ymin = 0,
  ymax = 1
)
rect_sampling_point_long <- c(40)
rect_sampling_point_long <- data.frame(
  xmin = rect_sampling_point_long,
  xmax = rect_sampling_point_long + 200,
  ymin = .2,
  ymax = .6
)

singe_time_point <- c(300)
singe_time_point <- data.frame(
  xmin = singe_time_point,
  xmax = singe_time_point + 10,
  ymin = 0,
  ymax = 1
)

ggplot() +
  geom_rect(data=rect_sampling_point_unique, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill='blue', alpha=0.2) +
  geom_rect(data=rect_sampling_point_long, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill='green', alpha=0.4) +
  geom_rect(data=singe_time_point, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill='black', alpha=0.4) +
  geom_line(data=smooth_yt, aes(y=V1, x=x, col='seroprevalence')) +
  ylab('seroprevalence') +
  xlab('julian days') +
  #ggtitle('sampling schemes of henipavirus/filovirus seroprevalence studies present in published literature') +
  labs(caption = "Sampling schemes present in published studies of henipavirus/filoviru sereprevalence in bats.
                  Black represents a single sampling event.
                  Blue represents a longitudinal sampling at bimonthly intervals. 
                  Green represents a potentially longitudinal sampling, but the study reports only a single mean estimate of seroprevalence for a given sample period.", parse = TRUE)


#.............ehhhh just gonna add on that alex shit rn and try and figure it out.........
plot(smooth_y[120:280], type='l')

y<-smooth_y[120:280]

smooth_yt <- as.data.frame((as.matrix(y))) %>%
  mutate(x = row_number()) %>%
  #mutate(x1 = (1.75*pi*(x/161))) %>%
  #mutate(x2=.3*sin(x1-1.3)+.3) %>%
  mutate(x1 = (2*pi*(x/161))) %>%
  mutate(x2=sin(x1)) %>%
  mutate(x3=x2)

plot(smooth_yt$x3, type='l')
lines(smooth_yt$V1, type='l')

smooth_yt %>%
  ggplot(aes(x=x3, y= V1))+
  geom_point() +
  geom_line()

lm(V1~x2,data=smooth_yt) %>% summary()

#.........eh draw a circle i guess...............................
library(plotrix)

x<-seq(0,365,by=1)
amp <- rnorm(n=length(x), mean=0, sd=.1)
y<-sin(.04*x)+amp
y<-y*sd(y)
plot(y, type='l')

y<-y[0:160]
plot(y, type='l')

smooth_yt <- as.data.frame((as.matrix(y))) %>%
  mutate(V2 = V1 + abs(min(smooth_yt$V1))) %>%
  mutate(x = row_number()) %>%
  mutate(x1 = (2*pi*(x/161))) %>%
  mutate(x2=sin(x1))

plot(smooth_yt$x2)
lines(smooth_yt$V1)
plot(x=smooth_yt$x2, y = smooth_yt$V1)

par(mar = c(4,4,4,4)) 
plot(x=0,y=0,xlim=c(-1.5,1.5), ylim=c(-1.5,1.5))
plot(x=0,y=0)
draw.circle(y=0,x=0,radius=mean(smooth_yt$V2))
points(x = smooth_yt$x2, y = smooth_yt$V1)


n<-1000
plot.new()
frame()
x<-runif(n,-1,1)
y<-runif(n,-1,1)
for (i in 1:n) { plot(x[i],y[i])}
draw.circle(0,0,1,nv=1000,border=NULL,col=NA,lty=1,lwd=1)





