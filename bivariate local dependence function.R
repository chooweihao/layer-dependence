#test "marginal correlation" for various copulas

#setwd("C:\Users\sgarcwh\Documents\PhD\marginal dependence")

#normal copula
norm <- normalCopula(0.75,dim=2)
u <- rCopula(100000, norm)

#gumbel copula
gum <- gumbelCopula(2,dim=2)
u <- rCopula(100000, gum)

#frank copula
fra <- frankCopula(10,dim=2)
u <- rCopula(100000, fra)

#joe copula
joe <- joeCopula(param=5, dim = 2)
u <- rCopula(100000, joe )


#t copula
t <- tCopula(0.75,df=200,dim=2)
u <- rCopula(100000, t)


#clayton copula
clay <- claytonCopula(3,dim=2)
u <- rCopula(100000, clay)


####################


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; rho=cor(u1,u2);
n=100
cutoff=matrix((1:n)/n);
marcor=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor[i]=cov(u2,ind1)/cov(u1,ind1);
}


#compute bivariate local dependence function

k=10000;
s=matrix(1:k)/k;
t=0; difft=0; adjdifft=0;

for (i in 1:k)
{
t[i]=log(dCopula(c(s[i],0.99),clay )/dCopula(c(s[i],0.01),clay ));
}

for (i in 1:k-1)
{
difft[i]=(t[i+1]-t[i])*k;
adjdifft[i]=2*exp(difft[i])/(1+exp(difft[i]))-1;
}

#plot(s[100:k-100],difft[100:k-100])
#plot(s[100:k-100],adjdifft[100:k-100])



#setEPS();
#postscript("biclayton.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u,alpha))),ylab=expression(paste(list(v,rho[alpha],gamma[alpha]))),cex.lab=2);
lines(s[1000:k-1000],adjdifft[1000:k-1000]^5000,lty=1,lwd=4,col="red");
lines(cutoff,marcor,lty=1,lwd=4,col="blue");
#dev.off();

