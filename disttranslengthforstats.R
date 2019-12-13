#using Sim.DiffProc and fitsde
rm(list=ls())
library("Sim.DiffProc")

#parameters
p=12
h=1.0
r=.6
b=1.0
sigma=0.05
start=1.13
#start=0.83

library(rootSolve)
# the model : lake

#create simulated data
f <- expression((r*x^p)/(x^p+h^p)+a-b*x)
g <- expression(sigma)

rate <- function(x, a = 0.67) (r*x^p)/(x^p+h^p)+a-b*x
# find all roots within the interval [0,10]
Eq   <- uniroot.all(rate, c(0, 2))
start=Eq[3]+rnorm(1,0,0.01)

ap=seq(0.64254962, 0.55, -0.01)
reps=100
freqT=1000
lenT=100
tcon=lenT/freqT
dval=0.01

tr = matrix(nrow=reps,ncol=length(ap))
eqd=vector()
lg=vector()


for(j in 1:length(ap)){
  plot(0, ylim=c(0.4,1.3), xlim=c(0,50))
  print(ap[j])
  a = ap[j]
rate <- function(x, a = ap[j]) (r*x^p)/(x^p+h^p)+a-b*x
EqT   <- uniroot.all(rate, c(0, 2))
ed=EqT[1]

eqd[j]=ed
print(ed)
TEX.sde(object=c(drift = f, diffusion = g))
f <- expression((r*x^p)/(x^p+h^p)+a-b*x)
for (i in 1:reps) {
HWV <- snssde1d(drift=f,diffusion=g,x0=start,N=freqT,T=lenT)
tr[i,j]=min(which((abs(HWV$X-ed)<dval)))
}
}

lt<-apply(tr, 2, function(x) sum(is.infinite(x)))
trn<-as.data.frame(ifelse(is.infinite(tr), freqT, tr))
mt<-apply(trn, 2, function(x) mean(x, na.rm=TRUE))
mt<-mt*tcon
mxt<-apply(trn, 2, function(x) max(x, na.rm=TRUE))
mxt<-mxt*tcon
mnt<-apply(trn, 2, function(x) min(x, na.rm=TRUE))
mnt<-mnt*tcon
trn<-rev(trn)*tcon
colnames(trn)<-rev(ap)


for (j in 1:length(ap)){
  lg[j]<-sum(tr[,j]>(20/tcon))/reps*100
}

boxplot((trn), names=rev(round(ap,2)), log = "y",xlab="bifurcation parameter: phosophorus loading (a)", ylab="length of transient (yrs) ",
        main=(paste("Scheffer's lake model with white noise (sigma=", sigma, "; ",reps, "replicates)")))
abline(v=length(ap)+.4, col="red", lwd=2)
text(length(ap)+.6, 50, "Bifurcation", srt=90)
library(latex2exp)
#text(3,60, TeX('$dX_{t}=\\left( \\frac{ r X_{t}^p}{(X_{t}^p + h^p)} - b X_{t} + a \\right)dt +\\sigma dW_{t}$'))
text(3,60, TeX('$dX_{t}=\\left( \\frac{ 0.6 X_{t}^{12}}{(X_{t}^{12} + 1.0^{12})} - 1.0 X_{t} + a \\right)dt +\\0.05 dW_{t}$'))

