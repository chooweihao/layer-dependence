#test "marginal correlation" for various copulas

#setwd("C:\Users\sgarcwh\Documents\PhD\marginal dependence")

#normal copula
norm <- normalCopula(0.5,dim=2)
u <- rCopula(10000000, norm)

#gumbel copula
gum <- gumbelCopula(5,dim=2)
u <- rCopula(10000000, gum)


#frank copula
fra <- frankCopula(5,dim=2)
u <- rCopula(10000000, fra)


#joe copula
joe <- joeCopula(param=5, dim = 2)
u <- rCopula(10000000, joe )


#t copula
t <- tCopula(0.75,df=200,dim=2)
u <- rCopula(10000000, t)


#clayton copula
clay <- claytonCopula(3,dim=2)
u <- rCopula(10000000, clay)

#structural copula
n=10000;
x=rgamma(n,1,1); e1=rnorm(n,0,1); e2=rnorm(n,0,1);
z1=x+e1; z2=x+e2;
u1=rank(z1)/n; u2=rank(z2)/n;
u=cbind(u1,u2);

#mixing copula
n=1000000; a=0.5; 
e1=runif(n,0,1); e2=runif(n,0,1); x=runif(n,0,1); z=runif(n,0,1);
u1=(z<=a)*x+(z>a)*e1;
u2=(z<=a)*x+(z>a)*e2;
u=cbind(u1,u2);

####################


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; rho=cor(u1,u2);
n=50
cutoff=matrix((0:n)/n);
marcor=0; #marginal correlation
conexp=0; #conditional mean
gconexp=0; #gradient
convar=0; #conditional variance
midpoint=0; #midpoint
corcurve=0; #correlation curve
uvar=var(u2);

for (i in 1:n)
{
k=cutoff[i+1]; j=cutoff[i]; ind=(u1<=k)*(u1>j);
conexp[i]=mean(u2[ind==1]); 
convar[i]=var(u2[ind==1]);
midpoint[i]=(k+j)/2;
}

for (i in 2:n)
{
gconexp[i]=(conexp[i]-conexp[i-1])/(midpoint[i]-midpoint[i-1]);
corcurve[i]=gconexp[i]/(gconexp[i]^2+convar[i]/uvar)^0.5;
}

for (i in 2:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor[i]=cov(u2,ind1)/cov(u1,ind1);
}

#setEPS();
#postscript("cfrank.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u,alpha))),ylab=expression(paste(list(v,c[alpha],rho[alpha]))),cex.lab=2);
lines(midpoint[2:n],corcurve[2:n],lty=1,lwd=4,col="red");
lines(cutoff[2:n],marcor[2:n],lty=1,lwd=4,col="blue");
#lines(midpoint,conexp,lty=2,lwd=4);
#lines(midpoint,convar,lty=2,lwd=4);
#lines(midpoint[2:n],gconexp[2:n],lty=2,lwd=4);
#dev.off();

