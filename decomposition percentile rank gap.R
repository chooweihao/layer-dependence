#test "marginal correlation" for various copulas

#setwd("C:\Users\sgarcwh\Documents\PhD\marginal dependence")

#normal copula
norm <- normalCopula(0.75,dim=2)
u <- rCopula(100000, norm)

#gumbel copula
gum <- gumbelCopula(2,dim=2)
u <- rCopula(100000, gum)
plot(u,cex=0.1);

#frank copula
fra <- frankCopula(10,dim=2)
u <- rCopula(100000, fra)
plot(u,cex=0.1);

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


#compute decompositions
data=u; u1=data[,1]; u2=data[,2]; rho=cor(u1,u2);
n=100
cutoff=matrix((1:n)/n);
prg=0; #percentile rank gap
conprob=0; #concordance probability
disp=0; #dispersion

for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
prg[i]=cov(u2,ind1)/cov(u1,ind1);
conprob[i]=(mean((u1<=k)*(u2<=k))-k^2)/(k*(1-k));
disp[i]=2*mean(abs(u1-u2)*((u1<=k)*(u2>k)+(u1>k)*(u2<=k)))/mean((u1<=k)*(u2>k)+(u1>k)*(u2<=k));
}

setEPS();
postscript("frankdecom.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u,alpha))),ylab=expression(paste(list(v,rho[alpha]))),cex.lab=2);
lines(cutoff,conprob,lty=1,lwd=4,col="red");
lines(cutoff,disp,lty=1,lwd=4,col="blue");
lines(cutoff,1-(1-conprob)*disp,lty=2,lwd=4);

dev.off();

