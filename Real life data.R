#market returns
mydata = read.csv("nasdaq ftse returns v2.csv");

nasdaq=mydata[,1];
ftse=mydata[,2];

n=length(nasdaq);

nasdaq1=rank(nasdaq)/n;
ftse1=rank(ftse)/n;

rhodata=cor(ftse1,nasdaq1);


u=cbind(nasdaq1,ftse1);
plot(u);



#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; rho=cor(u1,u2);
n=100
cutoff=((1:n)/n);
marcor=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); ind2=(u1>k);
marcor[i]=cov(u2,ind1)/cov(u1,ind1);
}


setEPS();
postscript("NASDAQvsFTSE.eps")
par(mar=c(5,5,1,1));
plot(data,type="p",cex=0.25,xlab="NASDAQ",ylab="FTSE",cex.lab=2);
lines(cutoff,marcor,lty=1,lwd=4,col="red");
lines(cutoff,rho*rep(1,n),lty=2,lwd=4);
dev.off();


rough_target=marcor;


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

n=50000;
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

adj=1;

for (i in 1:20)
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






#plot fitted copula
w=floor(50000*runif(7418,0,1))+1;


setEPS();
postscript("NASDAQvsFTSE_fitted.eps")
par(mar=c(5,5,1,1));
plot(us[w],vs[w],type="p",cex=0.25,xlab="NASDAQ",ylab="FTSE",cex.lab=2);
lines((1:99)/100,predict(target),lty=1,lwd=4);
lines((1:100)/100,rough_target,lty=1,lwd=4,col="red");
dev.off();


#gaussian fitted

#normal copula
norm <- normalCopula(0.415,dim=2)
u <- rCopula(100000, norm)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; rho=cor(u1,u2);
n=100
cutoff=((1:n)/n);
marcor=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); ind2=(u1>k);
marcor[i]=cov(u2,ind1)/cov(u1,ind1);
}


setEPS();
postscript("gaussianfitted.eps")
par(mar=c(5,5,1,1));
plot(u1[1:7418],u2[1:7418],type="p",cex=0.25,xlab="NASDAQ",ylab="FTSE",cex.lab=2);
lines(cutoff,marcor,lty=1,lwd=4);
lines((1:100)/100,rough_target,lty=1,lwd=4,col="red");
dev.off();

#plot densities

varf=rhodata/(1-rhodata);
z=(-200:500)/100;


setEPS();
postscript("densities.eps")
par(mar=c(5,5,1,1));
plot(density(Fz),col="red",xlab="s",ylab="Density of s",cex.lab=2,main="",lwd=4);
lines(z,dnorm(z,mean(Fz),varf^0.5),lwd=4);

dev.off();


