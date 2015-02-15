
####################



#rho=1
c=0.75; 
fra <- frankCopula(6,dim=2)
u <- rCopula(100000*c, fra)
u1a=c*u[,1]; u2a=c*u[,2];
fra <- frankCopula(6,dim=2)
u <- rCopula(100000*(1-c), fra)
u1b=c+(1-c)*u[,1]; u2b=c+(1-c)*u[,2];
u1=cbind(u1a,u1b);
u2=cbind(u2a,u2b);


i=floor(runif(50000,0,1)*40000)+1;
u1=u1[i];
u2=u2[i];

plot(u1,u2,cex=0.1);


#plot initial graph
#setEPS();
#postscript("perfect.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u))),ylab=expression(paste(list(v))),cex.lab=2);
lines(0.5*rep(1,100),(1:100)/100,lty=2,lwd=3);
lines((1:100)/100,0.5*rep(1,100),lty=2,lwd=3);
#dev.off();








#frank copula
fra <- frankCopula(10,dim=2)
u <- rCopula(1000000, fra)
u1=u[,1]; u2=u[,2];

#setEPS();
#postscript("imperf.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u))),ylab=expression(paste(list(v))),cex.lab=2);
lines(0.75*rep(1,100),(1:100)/100,lty=2,lwd=3);
lines((1:100)/100,0.75*rep(1,100),lty=2,lwd=3);
#dev.off();

a=0.75;
cov(u2,(u1>a))/cov(u1,(u1>a));


cor((u1>a),(u2<=a));
mean(abs(u1-u2)*((u1-a)*(u2-a)<0))/mean((u1-a)*(u2-a)<0);
####################



