#test "marginal correlation" for various copulas

#setwd("C:\Users\sgarcwh\Documents\PhD\marginal dependence")

#normal copula
norm <- normalCopula(0.75,dim=2)
u <- rCopula(100000, norm)

u1=u[,1];
u2=u[,2];


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
u <- rCopula(1000000, clay)


####################


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; rho=cor(u1,u2);
n=500
cutoff=matrix((1:n)/n);
marcor=0; #marginal correlation
marcor1=0; 
marcor2=0; 
marcor3=0;
c1=0.25; c2=0.5; c3=0.75;
n1=5; n2=10; n3=20;

for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor[i]=cov(u2,ind1)/cov(u1,ind1);
marcor1[i]=cov((u2>c1),ind1)/cov((u1>c1),ind1);
marcor2[i]=cov((u2>c2),ind1)/cov((u1>c2),ind1);
marcor3[i]=cov((u2>c3),ind1)/cov((u1>c3),ind1);
}

for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor[i]=cov(u2,ind1)/cov(u1,ind1);
marcor1[i]=cov(u2^n1,ind1)/cov(u1^n1,ind1);
marcor2[i]=cov(u2^n2,ind1)/cov(u1^n2,ind1);
marcor3[i]=cov(u2^n3,ind1)/cov(u1^n3,ind1);
}

setEPS();
postscript("eml.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u,alpha))),ylab=expression(paste(list(v,rho[alpha]^phi))),cex.lab=2);
lines(cutoff,marcor,lty=2,lwd=5,col="black");
lines(cutoff,marcor1,lty=1,lwd=4,col="blue");
lines(cutoff,marcor2,lty=1,lwd=4,col="red");
lines(cutoff,marcor3,lty=1,lwd=4,col="darkgreen");
legend(0.4,0.3,c(expression(n==5),expression(n==10),expression(n==20))
,lty=c(1,1),lwd=c(4,4),col=c("blue","red","darkgreen"), cex=2);
dev.off();


#transformed scale
n=n2;

setEPS();
postscript("eml2.eps")
par(mar=c(5,5.5,1,1));
plot(u1[1:1000]^n,u2[1:1000]^n,type="p",cex=1,
xlab=expression(paste(list(u^n,alpha^n)))
,ylab=expression(paste(list(v^n,rho[alpha],rho[alpha]^phi))),cex.lab=2);
lines(cutoff^n,marcor,lty=2,lwd=5,col="black");
#lines(cutoff^n,marcor1,lty=1,lwd=4,col="blue");
lines(cutoff^n,marcor2,lty=1,lwd=4,col="red");
#lines(cutoff^n,marcor3,lty=1,lwd=4,col="darkgreen");
dev.off();





