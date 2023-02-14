#Set up the dataset
x.prime <- 1/c(1.0,1.5,2.0,3.0,4.0,5.0,6.0,7.5,8.5,10.0,12.5,15.0,17.5,20.0,25.0,30.0,35.0,40.0)
y.prime <- 1/c(2.1,2.5,4.9,5.5,7.0,8.4,9.6,10.2,11.4,12.5,13.1,14.6,17.0,16.8,18.6,19.7,21.3,21.6)
x <- c(1.0,1.5,2.0,3.0,4.0,5.0,6.0,7.5,8.5,10.0,12.5,15.0,17.5,20.0,25.0,30.0,35.0,40.0)
y <- c(2.1,2.5,4.9,5.5,7.0,8.4,9.6,10.2,11.4,12.5,13.1,14.6,17.0,16.8,18.6,19.7,21.3,21.6)
df <- data.frame(x.prime, y.prime, x, y)

#1.a Fit a linear regression model and find starting values
lm1 <- lm(df$y.prime~df$x.prime)
lm1
g0 <- 1/coef(lm1)[1]
g1 <- coef(lm1)[2]/coef(lm1)[1]

#1.b Find the least square estimates
start <- nls(y ~ gamma0*x/(gamma1+x), 
                 data = df, 
                 start = list(gamma0=g0, gamma1=g1))
gamma.hat <- start$m$getAllPars()
gamma.hat

#2.a Plot the data and the fit
plot(df$x, df$y, xlab="concentration", ylab="velocity")
lines(df$x, start$m$fitted(), col="red")

#2.b Draw the residual plot and QQ plot
plot(start$m$fitted(), start$m$resid(), xlab="Fitted values", ylab="Residuals"); 
abline(h=0, col="red", lty=2)

qqnorm(start$m$resid(), ylab="Residuals")
qqline(start$m$resid(), col="red")

#2.c Test gamma1 = 20
MSE <- sum(start$m$resid()^2)/16
Var <- MSE*solve(t(J)%*%J)
Var
s.gamma1 <- sqrt(Var[2,2])
t <- (gamma.hat[2]-20)/s.gamma1
t

#3 Apply Bootstrapping and calculate the confidence intervals
n=dim(df)[1]

gamma0.v=NULL
gamma1.v=NULL

set.seed(1)

for(i in 1:1000)
{
  newID = sample(1:n,replace=T)
  new.df = df[newID,]
  newfit =  nls(y ~ gamma0*x/(gamma1+x), 
                data = new.df, start = list(gamma0=g0, gamma1=g1))
  gamma0.v=c(gamma0.v,coef(newfit)[1])
  gamma1.v=c(gamma1.v,coef(newfit)[2])
}
quantile(gamma0.v, probs = c(0.025, 0.975))
quantile(gamma1.v, probs = c(0.025, 0.975))

#Compare to the confidence interval from the Large Sample Theory
J <- start$m$gradient()
sigma2 <- sum(start$m$resid()^2)/(nrow(J)-ncol(J))
se.gamma0 <- sqrt(sigma2)*sqrt( solve(t(J)%*%J)[1,1])
se.gamma1 <- sqrt(sigma2)*sqrt( solve(t(J)%*%J)[2,2])
CI.gamma0 <- gamma.hat[1] + c(-1,1)*se.gamma0*qt(p = 0.975, df = nrow(J)-ncol(J))
CI.gamma1 <- gamma.hat[2] + c(-1,1)*se.gamma1*qt(p = 0.975, df = nrow(J)-ncol(J))

CI.gamma1
