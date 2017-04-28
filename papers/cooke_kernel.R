library(plot3D)
library(np)
DAT.NP <- read.csv("~/R/DAT.csv")
K <- DAT.NP$K
L <- DAT.NP$L
Y <- DAT.NP$Y
n <- 1000
KdY <- K/Y
LdY <- L/Y
h1<-1.06*sd(K)*n^(1/6) 
h2<-1.06*sd(L)*n^(1/6) 
## Question 1a
bw1 <- npregbw(formula=Y~K+L,regtype="lc",bandwidth.compute=FALSE,bws=c(h1,h2))
lcls <- npreg(bws = bw1, gradients = TRUE)
summary(lcls)
plot (bw1)

