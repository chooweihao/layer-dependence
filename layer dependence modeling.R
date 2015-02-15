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




#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; rho=cor(u1,u2);
n=100
cutoff=((1:n)/n);
marcor=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); ind2=(u1>k);
marcor[i]=cov(u2,ind1)/cov(u1,ind1);
#2*(mean(u2*ind2)/mean(ind2)-mean(u2*ind1)/mean(ind1));
}

#setEPS();
#postscript("normal.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u,alpha))),ylab=expression(paste(list(v,rho[alpha]))),cex.lab=2);
lines(cutoff,marcor,lty=1,lwd=5,col="red");
#lines(cutoff,rho*rep(1,n),lty=2,lwd=5);
#dev.off();





x=cutoff;
y=marcor;
target=nls(y~a+b*x+c*x^2+d*x^3+e*x^4,start=list(a=0,b=0,c=0,d=0,e=0));
summary(target);





#function to calculate layer dependence


layerdep=function(u,v)
{
n=100; a=c((1:n)/n); rhoa=0; 
for (i in 1:n){ind=(u<=a[i]); rhoa[i]=cov(v,ind)/cov(u,ind);}
return(rhoa);
}

layerdep1=function(u,v)
{
n=100; a=c((1:n)/n); rhoa=0; 
for (i in 1:n){ind=(u<=a[i]); rhoa[i]=cov(v,ind)/cov(u,ind);}
plot(u,v,type="p",cex=0.1);
lines(a,rhoa,lty=1,lwd=5,col="red");
}



#iteration

n=20000;
z=c(1:(n-1))/n;
Fz=rho*qnorm(z,0,1);
error=rnorm(n-1,0,1);
error1=rnorm(n-1,0,1);

us=z;
vs=rank(Fz+error)/(n-1);
plot(us,vs,cex=0.1);
rhoxx=layerdep(us,vs);
xx=(1:100)/100;
fitted=nls(rhoxx~a+b*xx+c*xx^2+d*xx^3+e*xx^4,
start=list(a=0,b=0,c=0,d=0,e=0));


hist(Fz,25);

adj=2;

for (i in 1:10)
{
Dz=diff(Fz,lag=1,differences=1);

for (j in 2:(n-1))
{
t=predict(target,newdata=list(x=j/n));
f=predict(fitted,newdata=list(xx=j/n));
Fz[j]=Fz[j-1]+Dz[j-1]*(t/f)^adj;

}
hist(Fz,25);
us=rank(Fz+error)/(n-1);
vs=rank(Fz+error1)/(n-1);
rhoxx=layerdep(us,vs);
fitted=nls(rhoxx~a+b*xx+c*xx^2+d*xx^3+e*xx^4,
start=list(a=0,b=0,c=0,d=0,e=0));
}


layerdep1(us,vs);
lines((1:99)/100,predict(target),lty=2,lwd=5,col="green");


