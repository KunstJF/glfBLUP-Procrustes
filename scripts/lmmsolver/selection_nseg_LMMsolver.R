# Martin Boer, 24 juli 2024

library(LMMsolver)

set.seed(1234)

days <- c(0, 11, 19, 25,28 ,32 ,36 ,39 ,48 ,50 ,53, 63,
          68, 74, 80 , 88 ,95, 102, 109, 118, 120)

x0 <- seq(0,120,by=0.1)

f <- function(x) { exp(-0.01*(x-60)^2) + exp(-0.01*(x-20)^2) +
                                         exp(-0.01*(x-100)^2) }
#f <- function(x) { sin(2*pi*x/20) }

nrepl <- 25
d <- rep(days,each=nrepl)

plot(x0,exp(-0.001*(x0-60)^2))

y0 <- f(x0)
y_true <- f(d)
y <- y_true + rnorm(n=21*nrepl,sd=0.1)
dat <- data.frame(y=y,x=d,y_true=y_true)

# 1. number of segments nseg=7 too low, overfitting:
obj <- LMMsolve(y~1,spline=~spl1D(x, nseg=7), data=dat)
summary(obj)
pred <- obtainSmoothTrend(obj,grid=500,includeIntercept = TRUE)
plot(x=dat$x,y=dat$y,col='red',xlab='time',ylab='y')
lines(x=pred$x,y=pred$ypred,lw=2.5,col='green')
lines(x=x0,y=y0,col='blue')
abline(h=0)

# 2. logL as function of number of segments
nseg <- c(7:50)
logL <- rep(NA, length(nseg))
for (i in 1:length(nseg)) {
  obj <- LMMsolve(y~1,spline=~spl1D(x,nseg=nseg[i]), data=dat)
  logL[i] <- obj$logL
}
plot(nseg,y=logL,type='b')

# 3. extra test for non-smooth effect with many segments:
# in this example the non-smooth effect xF goes to zero, as expected for
# the simulated smooth data
dat$xF <- as.factor(dat$x)
obj <- LMMsolve(y~1,spline=~spl1D(x, nseg=50), random=~xF, data=dat)
summary(obj)


