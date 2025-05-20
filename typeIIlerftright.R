rm(list=ls())
set.seed(NULL)
theta=1
ds=10000
n=30;r1=6;r2=1
cdf=function(par,x)
{
  a=exp(-(par^2/x^2))
  return(a)
}

QAD=function(par){
  if (r1<r2) rr=r2 else rr=r1
  j=c((rr+1):(n-rr))
  # x=sort(xs)
  xad=xs[j]
  adx=sort(xad,decreasing=TRUE)
  AD=(-n-(1/n)*sum((2*j-1)*(log(cdf(par,xad))+log(1-cdf(par,adx)))))#n-i+1
  return(AD)
}

QCV=function(par){
  j=c((r1+1):(n-r2))
  #x=sort(xs)
  CV=1/(12*n)+sum((cdf(par,x)-(j-0.5)/n)^2) #yeni
  return(CV)
}

QEKK=function(par){
  j=c((r1+1):(n-r2))
  #x=sort(xs)
  Q=3#yeni
  return(Q)
}
#
####Pdf of IRD####

pdfIR=function(x,theta){(2*theta^2/x^3)*exp(-(theta^2/x^2))}

####Cdf of IRD####
cdfIR=function(x,theta){exp(-(theta^2/x^2))}
gz1=function(z){z*pdfIR(z,1)/(1-cdfIR(z,1))}
gz2=function(z,the){z*pdfIR(z,the)/(1-cdfIR(z,the))}

####Likelihood function####

lf=function(theta)
{
  res=r1*(log(cdfIR(x[1],theta)))+
    r2*log(1-(cdfIR(x[xl],theta)))+
    sum(log(pdfIR(x[(1):(xl)],theta)))
  return(-res)
}
sonuc=NULL;mse=NULL
for(jj in 1:ds)
{
  cat("\14",jj)
  u=runif(n)
  y=1/(-log(u))^(1/2)*theta
  xs=sort(y)
  x=xs[(r1+1):(n-r2)]
  xl=length(x)
####MMLE1####
  ln=function(x)log(x)
  h=1/((-log((n-r2)/(n+1)))^(1/2))
  K2=-h^2*(2*exp(-1/(h^2))*(-3*h^2+3*h^2*exp(-1/(h^2))+2)/h^6/(-1+exp(-1/(h^2)))^2)
  # KS=-2/h^4*exp(-1/(h^2))*(-3*h^2+3*h^2*exp(-1/(h^2))+2)/(-1+exp(-1/(h^2)))^2
  A=sum(1/x^2)
  za=x[1]
  xa=za
  xb=max(x)
  A0=2*(n-r1-r2)
 #thetammle3=1/2*2^(1/2)*((A*xa+r1)*xa*(A0+r2*K2))^(1/2)/(A*xa+r1)#eski hali
  mmle1=(xa^2*(A0+r2*K2))^(1/2)*1/sqrt(2*(A*xa^2+r1))#thetammle3 sade hali düzeltilmiş
  lastmmle=1/2*2^(1/2)/(A*xa^2+r1)^(1/2)*(A0+r2*K2)^(1/2)*xa #benim bulduğum 
  ####MMLE2####
  zz = 1/((-log((n-r2)/(n+1)))^(1/2))
  
  mmle2=1/2*2^(1/2)/(A*xa^2+r1)^(1/2)*((A0+r2*gz1(zz)))^(1/2)*xa
  #st=1/(2*(A*xa^2+r1))^(1/2)*((A0+r2*gz1(zz))*xa^2)^(1/2)#manuscripte yazdığım şekli
  thet0=mmle2
  ####UPDATE MMLE2####
  for(i in 1:10){
    thet2=1/2*2^(1/2)/(A*xa^2+r1)^(1/2)*(A0+r2*gz2(xb,thet0))^(1/2)*xa
    thet0=thet2
  }
  ####MLE, AD, CVM, EKK####
  mle=try(optim(c(theta),lf),silent = T)$par
  ad=try(optim(c(theta),QAD),silent = T)$par
  cvm=try(optim(c(theta),QCV),silent = T)$par
  ekk=try(optim(c(theta),QEKK),silent = T)$par
  #wekk=try(optim(c(theta),QWEKK),silent = T)$par
  
  
  sonuc=rbind(sonuc,c("mle"=mle,"mmle1"=mmle1,"mmle2"=mmle2,"ad"=ad,"cvm"=cvm,"ekk"=ekk))#mmle1 Klı olan
}
bias=apply(sonuc,2, mean)-theta
mse_mle=mean((sonuc[,1]-theta)^2) #mle
mse_mmle1=mean((sonuc[,2]-theta)^2) #K olan mmle
mse_mmle2=mean((sonuc[,3]-theta)^2) #bizim mmle
mse_ad=mean((sonuc[,4]-theta)^2)
mse_cvm=mean((sonuc[,5]-theta)^2)
mse_ekk=mean((sonuc[,6]-theta)^2)
mse=c("mle"=mse_mle,"mmle1"=mse_mmle1,"mmle2"=mse_mmle2,"ad"=mse_ad,"cvm"=mse_cvm,"ekk"=mse_ekk)
mse_all=round(c(mse),4);mse_all
bias_all=round(c(bias),4);bias_all

res=c("theta"=theta,"n"=n,"r1"=r1,"r2"=r2,"mse"=mse_all,"bias"=bias_all);res

 library(xlsx)
write.xlsx(res, file="30_6_1.xlsx")
# c(theta,mle,mmle1,lastmmle,mmle2)
