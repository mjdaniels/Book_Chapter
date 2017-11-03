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
agegrid<-seq(5,18,length=50)
data=list(y=y,age=age,n=n,T=A,agegrid=agegrid)
#   
# inits<-function(){
#     list(
#       tau.y = 1,
#       tau.theta02 = 1,
#       xi0 = c(1,1),
#       nu0 = c(1,1),
#       omega0 = c(1,1),
#       alpha = 2
#     )	
# } 
##### parameters  

para<-c("A", "mu0", "tau.mu","sigmasq0.i","nu","mu.fit")
para<-c("A", "mu0", "tau.mu","sigmasq0.i","nu")
nlmefits3<-jags(data=data, inits=NULL, parameters.to.save=para, model.file="D:\\dandan\\short course\\examples\\dmd.dir\\nlme1_new.txt",
              n.chains=1, n.iter=10000, n.burnin=5000,n.thin=5)

print(nlmefits3)
nlme_results<-nlmefits3$BUGSoutput$sims.list
fitted<-apply(nlme_results,c(2,3),mean)
observed<-y[c(1,16,3,7,2,4,8,9),]
age_obs<-age[c(1,16,3,7,2,4,8,9),]
range(age_obs)


########### plot goodness of fit
par(mfrow=c(2,4))
plot(agegrid,fitted[1,],lty=1,type="l",ylab="Measure",xlab="age",xaxt='n',ylim=c(0,1),xlim=c(5,18))
axis(1,at=seq(5,18,by=1),labels=seq(5,18,by=1))
points(age_obs[1,],observed[1,])
for (i in 2:8) {
  plot(agegrid,fitted[i,],lty=1,type="l",ylab="Measure",xlab="age",xaxt='n',ylim=c(0,1),xlim=c(5,18))
  axis(1,at=seq(5,18,by=1),labels=seq(5,18,by=1))
  points(age_obs[i,],observed[i,])
}

lines(agegrid,colMeans(nlmefits3$BUGSoutput$sims.list$mupop),col="red",lty=1,type="l",lwd=2)
legend(5,1, c("subject-level","population-level    "), lty=c(1,1),lwd=c(1,2),col=c("black","red"))






save(nlme_results,file="D:\\dandan\\short course\\examples\\dmd.dir\\nlme_results.Rdata")
colMeans(cbind(nlme_results[[1]],nlme_results[[5]],nlme_results[[7]],
               nlme_results[[8]],nlme_results[[9]]))
1/colMeans(cbind(nlme_results[[1]],nlme_results[[5]],nlme_results[[7]],
                 nlme_results[[8]],nlme_results[[9]]))
estimate<-rbind(paste(round(mean(fits3$BUGSoutput$sims.list$A),2),"(",
                  round(quantile(fits3$BUGSoutput$sims.list$A,0.025),2),",",
                  round(quantile(fits3$BUGSoutput$sims.list$A,0.975),2),")",sep=""),
                paste(round(mean(1/fits3$BUGSoutput$sims.list$tau.y),1),"(",
                      round(quantile(1/fits3$BUGSoutput$sims.list$tau.y,0.025),1),",",
                      round(quantile(1/fits3$BUGSoutput$sims.list$tau.y,0.975),1),")",sep=""),
                paste(round(mean(fits3$BUGSoutput$sims.list$mu0),1),"(",
                  round(quantile(fits3$BUGSoutput$sims.list$mu0,0.025),1),",",
                  round(quantile(fits3$BUGSoutput$sims.list$mu0,0.975),1),")",sep=""),
            paste(round(mean(1/fits3$BUGSoutput$sims.list$tau.mu),1),"(",
                  round(quantile(1/fits3$BUGSoutput$sims.list$tau.mu,0.025),1),",",
                  round(quantile(1/fits3$BUGSoutput$sims.list$tau.mu,0.975),1),")",sep=""),
            paste(round(mean(fits3$BUGSoutput$sims.list$sigma0),1),"(",
                  round(quantile(fits3$BUGSoutput$sims.list$sigma0,0.025),1),",",
                  round(quantile(fits3$BUGSoutput$sims.list$sigma0,0.975),1),")",sep=""),
            paste(round(mean(1/fits3$BUGSoutput$sims.list$tau.sigma),1),"(",
                  round(quantile(1/fits3$BUGSoutput$sims.list$tau.sigma,0.025),1),",",
                  round(quantile(1/fits3$BUGSoutput$sims.list$tau.sigma,0.975),1),")",sep=""))
kable(cbind(estimate,estimate),format="latex")    



# mu<-NULL
# for (i in 1:133)
# {mu<-rbind(mu,apply(fits3$BUGSoutput$sims.list$mu.y[,i,],2,mean))}
# 
# par(mfrow=c(1,1))
# plot(time,mu[1,],lty=1,pch=1,type="o",ylab="Fat fraction in Sol",xlab="months",xaxt='n',ylim=c(0,0.7))
# axis(1,at=time,labels=time)
# for (i in 2:133) {
# lines(time,mu[i,],lty=i,pch=1,type="o")}



mu<-NULL
for (i in 1:10)
{mu<-rbind(mu,apply(nlmefits3$BUGSoutput$sims.list$mu.ypred[,i,],2,mean))}

####### plot individual vs population trajectories
par(mfrow=c(1,1))
plot(agegrid,mu[1,],lty=1,type="l",ylab="Measure",xlab="age",xaxt='n',ylim=c(0,1))
axis(1,at=seq(5,18,by=1),labels=seq(5,18,by=1))
for (i in 2:10) {
  lines(agegrid,mu[i,],lty=i,type="l")}

lines(agegrid,colMeans(nlmefits3$BUGSoutput$sims.list$mupop),col="red",lty=1,type="l",lwd=2)
legend(5,1, c("subject-level","population-level    "), lty=c(1,1),lwd=c(1,2),col=c("black","red"))