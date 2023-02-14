library(ROCR)
library(ResourceSelection)

x <- c(2,5,10,20,25,30)
n <- rep(500,6)
y <- c(72,103,170,296,406,449)
Y <- c(rep(0,500-72), rep(1,72), rep(0,500-103), rep(1,103), rep(0,500-170), rep(1,170), rep(0,500-296), rep(1,296),rep(0,94), rep(1,406), rep(0,51), rep(1,449))
X <- c(rep(2,500), rep(5,500), rep(10,500), rep(20,500), rep(25,500), rep(30,500))
n<- rep(1,3000)
data <- data.frame(cbind(x,y,n,prop))
dat <- cbind(X,Y,n)
prop <- y/n
dat <- data.frame(dat)

lm0 <- glm(prop ~ x, family = binomial(link = "logit"), data = data, weights = n)
summary(lm0)
plot(x,y/n, xlab = "size of deposit", ylab = "the probability a bottle will be returned")
lm1 <- glm(Y ~ X, family = binomial(link = "logit"), data = dat)
summary(lm1)
plot(data$x, data$prop, xlab="size of deposit", ylab="the probability a bottle will be returned")
lines(data$x, exp(data$x*lm1$coef[2]+lm1$coef[1])/(1+exp(data$x*lm1$coef[2]+lm1$coef[1])), col="red")
lines(data$x, lm3$fitted, col = "blue")

lm2 <- glm(Y~as.factor(X),family = binomial(link = "logit"), data = dat)
summary(lm2)

predicted = exp(lm1$coef[1]+lm1$coef[2]*x)/(exp(lm1$coef[1]+lm1$coef[2]*x)+1)
Expected <- cbind(500*(1-predicted),500*predicted)
Observed <- cbind(500-y, y)

Pearson.X2 <- sum((Observed-Expected)^2/Expected)
pchisq(q = Pearson.X2, df = 4, lower.tail = F)

Dev.X2 <- lm1$deviance - lm2$deviance
pchisq(Dev.X2, df = 4, lower.tail = F) 

x = c(0,10)
FI <- summary(lm0)$cov.scaled
se.x <-  sqrt(sum(x*(FI%*%x)))
beta.hat = lm0$coefficients
CI.logit <- sum(x*beta.hat) + se.x*c(qnorm(alpha/2), qnorm(1-alpha/2))
exp(CI.logit)

lm3 <- loess(prop~x, data = data, weights = n)
hoslem.test(dat$Y,fitted(lm1),g=6)

p <- predict(lm1, type="response")
pr <- prediction(p, dat$Y)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf, colorize = TRUE)
auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}
print(opt.cut(prf, pr))


