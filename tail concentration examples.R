#test "marginal correlation" for various copulas

#setwd("C:\Users\sgarcwh\Documents\PhD\marginal dependence")

#normal copula
norm <- normalCopula(0.75,dim=2)
u <- rCopula(200000, norm)

#gumbel copula
gum <- gumbelCopula(2,dim=2)
u <- rCopula(200000, gum)


#frank copula
fra <- frankCopula(10,dim=2)
u <- rCopula(200000, fra)


#joe copula
joe <- joeCopula(param=5, dim = 2)
u <- rCopula(100000, joe )


#t copula
t <- tCopula(0.75,df=200,dim=2)
u <- rCopula(100000, t)


#clayton copula
clay <- claytonCopula(3,dim=2)
u <- rCopula(50000, clay)
plot(u,cex=0.01);


####################


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; rho=cor(u1,u2);
n=100
cutoff=matrix((1:n)/n);
stdprob=0; #std prob
tailcon=0; #tail concentration
marcor=0; #percentile rank gap


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
stdprob[i]=(mean((u1<=k)*(u2<=k))-k^2)/(k*(1-k));
tailcon[i]=(k<=0.5)*mean((u1<=k)*(u2<=k))/mean(u1<=k)+(k>0.5)*mean((u1>k)*(u2>k))/mean(u1>k);
marcor[i]=cov(u2,ind1)/cov(u1,ind1);
}

setEPS();
postscript("vfrank.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u,alpha))),ylab=expression(paste(list(v,tau[alpha],rho[alpha]))),cex.lab=2);
lines(cutoff,tailcon,lty=1,lwd=4,col="blue");
lines(cutoff,stdprob,lty=2,lwd=5);
lines(cutoff,marcor,lty=1,lwd=4,col="red");
dev.off();

