#test "marginal correlation" for various copulas

#setwd("C:\Users\sgarcwh\Documents\PhD\marginal dependence")

#normal copula
norm <- normalCopula(0.5,dim=2)
u <- rCopula(100000, norm)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor1=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor1[i]=cov(u2,ind1)/cov(u1,ind1);
}


norm <- normalCopula(0.75,dim=2)
u <- rCopula(100000, norm)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor2=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor2[i]=cov(u2,ind1)/cov(u1,ind1);
}


norm <- normalCopula(0.9,dim=2)
u <- rCopula(100000, norm)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor3=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor3[i]=cov(u2,ind1)/cov(u1,ind1);
}






setEPS();
postscript("gaussianmul.eps")
par(mar=c(5,5,1,1));
plot(cutoff,marcor1,type="l",xlab=expression(alpha),ylab=expression(rho[alpha])
,lwd=4,ylim=c(0,1),col="red",cex.lab=2);
lines(cutoff,marcor2,lty=1,lwd=4,col="blue");
lines(cutoff,marcor3,lty=1,lwd=4,col="green");
legend(0.6,0.4,c(expression(rho==0.5),expression(rho==0.75),expression(rho==0.9))
,lty=c(1,1,1),lwd=c(4,4,4),col=c("red","blue","green"),cex=1.5);
dev.off();



#gumbel

gum <- gumbelCopula(2,dim=2)
u <- rCopula(100000, gum)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor1=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor1[i]=cov(u2,ind1)/cov(u1,ind1);
}


gum <- gumbelCopula(3,dim=2)
u <- rCopula(100000, gum)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor2=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor2[i]=cov(u2,ind1)/cov(u1,ind1);
}


gum <- gumbelCopula(4,dim=2)
u <- rCopula(100000, gum)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor3=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor3[i]=cov(u2,ind1)/cov(u1,ind1);
}






setEPS();
postscript("gumbelmul.eps")
par(mar=c(5,5,1,1));
plot(cutoff,marcor1,type="l",xlab=expression(alpha),ylab=expression(rho[alpha])
,lwd=4,ylim=c(0,1),col="red",cex.lab=2);
lines(cutoff,marcor2,lty=1,lwd=4,col="blue");
lines(cutoff,marcor3,lty=1,lwd=4,col="green");
legend(0.6,0.4,c(expression(theta==2),expression(theta==3),expression(theta==4))
,lty=c(1,1,1),lwd=c(4,4,4),col=c("red","blue","green"),cex=1.5);
dev.off();







#clayton

clay <- claytonCopula(2,dim=2)
u <- rCopula(100000, clay)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor1=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor1[i]=cov(u2,ind1)/cov(u1,ind1);
}


clay <- claytonCopula(3,dim=2)
u <- rCopula(100000, clay)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor2=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor2[i]=cov(u2,ind1)/cov(u1,ind1);
}


clay <- claytonCopula(4,dim=2)
u <- rCopula(100000, clay)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor3=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor3[i]=cov(u2,ind1)/cov(u1,ind1);
}






setEPS();
postscript("claytonmul.eps")
par(mar=c(5,5,1,1));
plot(cutoff,marcor1,type="l",xlab=expression(alpha),ylab=expression(rho[alpha])
,lwd=4,ylim=c(0,1),col="red",cex.lab=2);
lines(cutoff,marcor2,lty=1,lwd=4,col="blue");
lines(cutoff,marcor3,lty=1,lwd=4,col="green");
legend(0.6,0.4,c(expression(theta==2),expression(theta==3),expression(theta==4))
,lty=c(1,1,1),lwd=c(4,4,4),col=c("red","blue","green"),cex=1.5);
dev.off();




#frank

fra <- frankCopula(6,dim=2)
u <- rCopula(100000, fra)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor1=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor1[i]=cov(u2,ind1)/cov(u1,ind1);
}


fra <- frankCopula(8,dim=2)
u <- rCopula(100000, fra)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor2=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor2[i]=cov(u2,ind1)/cov(u1,ind1);
}


fra <- frankCopula(10,dim=2)
u <- rCopula(100000, fra)


#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; 
n=100
cutoff=matrix((1:n)/n);
marcor3=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor3[i]=cov(u2,ind1)/cov(u1,ind1);
}






setEPS();
postscript("frankmul.eps")
par(mar=c(5,5,1,1));
plot(cutoff,marcor1,type="l",xlab=expression(alpha),ylab=expression(rho[alpha])
,lwd=4,ylim=c(0,1),col="red",cex.lab=2);
lines(cutoff,marcor2,lty=1,lwd=4,col="blue");
lines(cutoff,marcor3,lty=1,lwd=4,col="green");
legend(0.6,0.4,c(expression(theta==6),expression(theta==8),expression(theta==10))
,lty=c(1,1,1),lwd=c(4,4,4),col=c("red","blue","green"),cex=1.5);
dev.off();






