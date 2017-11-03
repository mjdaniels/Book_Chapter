unlink(".RData")
rm(list=ls())

library(R2jags)
library(coda)
library(lattice)

################# read dataset
data<-read.table("D:\\dandan\\short course\\examples\\dmd.dir\\nlme.csv",sep=",",header=T)
y1<-data[order(as.numeric(data$id),data$time),]

id<-sort(unique(data$id))
time<-sort(unique(data$time))

n<-length(id)
A<-length(time)

y2<-data.frame(id[rep(1:n,each=A)],time[rep(1:A,n)])
colnames(y2)<-c("id","time")

y3<-merge(y2,y1,by=c("id","time"),all.x=T)
for (i in 2:(dim(y3)[1])){
	if (is.na(y3$age[i])==T) {
		y3$age[i]<-y3$age[i-1]+(y3$time[i]-y3$time[i-1])/12
	}
}

y4<-y3[order(y3$time,y3$id),]
y<-array(y4$sol,dim=c(n,A))
age<-array(y4$age,dim=c(n,A))

######################### Call jags run MCMC
K<-30
data=list(y=y,age=age,n=n,T=A,K=K)
  
##### initials
inits<-function(){
    list(
      tau.y = 1,
      mu00=15,
      tau.mu00=0.1,
      nu0=1,
      sigmasq00.i=0.05,
      A = 0.75,
      alpha = 2
    )	
} 
##### parameters  
para<-c("mu.y", "sigmasq.y", "A", "mu", "sigmasq", "theta", "xi.inv","nu","omega","theta0","xi0","nu0","omega0","alpha","cluster","group")
para<-c("A", "mu00", "tau.mu0","sigmasq0.i","nu")
nlmefits4<-jags(data=data, inits=NULL, parameters.to.save=para, model.file="D:\\dandan\\short course\\examples\\dmd.dir\\nlmeA1_new.txt",
              n.chains=1, n.iter=4000, n.burnin=3000,n.thin=1)
result<-nlmefits4$BUGSoutput$sims.matrix
print(nlmefits4)

