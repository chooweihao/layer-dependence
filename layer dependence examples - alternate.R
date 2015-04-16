#test "marginal correlation" for various copulas

#setwd("C:\Users\sgarcwh\Documents\PhD\marginal dependence")

#normal copula
norm <- normalCopula(0.617,dim=2)
u <- rCopula(100000, norm)

#gumbel copula
gum <- gumbelCopula(1.765,dim=2)
u <- rCopula(100000, gum)

#frank copula
fra <- frankCopula(4.45,dim=2)
u <- rCopula(100000, fra)

#joe copula
joe <- joeCopula(param=2.4, dim = 2)
u <- rCopula(100000, joe )

#t copula
t <- tCopula(0.68,df=1,dim=2)
u <- rCopula(100000, t)

#clayton copula
clay <- claytonCopula(1.51,dim=2)
u <- rCopula(100000, clay)

#structural copulas 1
n=100000; s=3.025;
e1=rnorm(n,0,1); e2=rnorm(n,0,1);
f=rgamma(n,1/2,1);
x=s*f+e1; y=s*f+e2;
u1=rank(x)/n; u2=rank(y)/n;
u=cbind(u1,u2);

#structural copulas 2
n=100000; s1=1.05; s2=1.05;
e1=rnorm(n,0,1); e2=rnorm(n,0,1);
f1=rgamma(n,1,1);
f2=rgamma(n,1,1);
x=s1*f1-s2*f2+e1; y=s1*f1-s2*f2+e2;
u1=rank(x)/n; u2=rank(y)/n;
u=cbind(u1,u2);

#structural copulas 3
n=50000; 
u1=0; u2=0;
c=rbeta(n,2,2);

for (i in 1:n)
{
z=runif(1,0,1);
u1[i]=(z<=c[i])*runif(1,0,c[i]) + (z>c[i])*runif(1,c[i],1);
u2[i]=(z<=c[i])*runif(1,0,c[i]) + (z>c[i])*runif(1,c[i],1);
}
u=cbind(u1,u2);
rho=cor(u1,u2);

#structural copulas 4
n=30000; 
u1=0; u2=0;

for (i in 1:n)
{
z=runif(1,0,1); s=runif(1,0,1); d=0.57;
c=(s<=d)*runif(1,0,1)+(s>d)*0.5;
u1[i]=(z<=c)*runif(1,0,c) + (z>c)*runif(1,c,1);
u2[i]=(z<=c)*runif(1,0,c) + (z>c)*runif(1,c,1);
}
u=cbind(u1,u2);
rho=cor(u1,u2);
plot(u,cex=0.01);

####################




#factor1 - gaussian
n=100000; 
e1=rnorm(n,0,1); e2=rnorm(n,0,1);
f=rnorm(n,0,1);
x=f+e1; y=f+e2;
u1=rank(x)/n; u2=rank(y)/n;
u=cbind(u1,u2);

#factor2 - chi square 1
n=100000; 
e1=rnorm(n,0,1); e2=rnorm(n,0,1);
f=rgamma(n,0.5,0.5);
x=f+e1; y=f+e2;
u1=rank(x)/n; u2=rank(y)/n;
u=cbind(u1,u2);



#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
u1=u1[2:99999]; u2=u2[2:99999];
n=100; cutoff=(1:n)/n;
rho=cor(u1,u2);
rho1=(2/3)^0.5*cor(log(u1/(1-u1)),u2);
rho2=(3^0.5)*cor(-log(1-u1),u2)-0.5*rho;
rho3=(3^0.5)*cor(log(u1),u2)-0.5*rho;



setEPS();
postscript("structural11.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u))),ylab=expression(paste(list(v))),cex.lab=2);
lines(cutoff,rho*rep(1,n),lty=1,lwd=5);
lines(cutoff,rho1*rep(1,n),lty=1,lwd=5,col="green");
lines(cutoff,rho2*rep(1,n),lty=1,lwd=5,col="blue");
lines(cutoff,rho3*rep(1,n),lty=1,lwd=5,col="red");
dev.off();










