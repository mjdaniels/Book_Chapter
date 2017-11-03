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

id.1<-matrix(rep(1:nrow(y),7),ncol=7)
y.1<-c(y)
age.1<-c(age)

y.2<-y.1[is.na(y.1)==FALSE]
age.2<-age.1[is.na(y.1)==FALSE]
id.2<-id.1[is.na(y.1)==FALSE]
nrec<-length(y.2)
nsub<-nrow(y)

grp<-matrix(0,nrec,nsub)
for (i in 1:nsub) {
  grp[id.2==i,i]<-1
}

data=list(y=y.2,age=age.2,id=id.2,nrec=nrec,
          nsub=nsub,grp=grp)

##### parameters  
para<-c("chi.diff","chi.diff1","logli")
nlmechk1<-jags(data=data, inits=NULL, parameters.to.save=para, model.file="D:\\dandan\\short course\\examples\\dmd.dir\\nlme1_chk.txt",
              n.chains=1, n.iter=20000, n.burnin=10000,n.thin=1)


####### model check
hist(nlmechk1$BUGSoutput$sims.list$chi.diff,xlab=expression(paste("difference in ",chi^2," discrepancy",sep="")),main="",breaks=seq(-200,200,by=10))
mean(nlmechk1$BUGSoutput$sims.list$chi.diff1)

##WAIC
logli<-as.matrix(nlmechk1$BUGSoutput$sims.list$logli)

lpd<-sum(log(colMeans(exp(logli))))
pwaic<-sum((apply(logli,2,sd))^2)
elpd<-lpd-pwaic
waic<- -2*elpd



################## cpo
logli1<-NULL

for (i in 1:133) {
  y.3<-y.2[id.2!=i]
  age.3<-age.2[id.2!=i]
  id.4<-id.2[id.2!=i]
  id.3<-id.4
  id.3[id.4>i]<-id.4[id.4>i]-1
  nrec.3<-length(y.3)
  nsub.3<-nsub-1
  ynew<-y.2[id.2==i]
  agenew<-age.2[id.2==i]
  nobs<-length(ynew)
  
  data=list(y=y.3,age=age.3,id=id.3,nrec=nrec.3,
            nsub=nsub.3,ynew=ynew,nobs=nobs,agenew=agenew)
  
  ##### parameters  
  para<-c("ynewlogli")
  #para<-c("munew","sigmasqnew.i","nu","sigmasq0.i")
  nlmecpo<-jags(data=data, inits=NULL, parameters.to.save=para, model.file="D:\\dandan\\short course\\examples\\dmd.dir\\nlme1_cpo.txt",
                n.chains=1, n.iter=4000, n.burnin=3000,n.thin=1)
  
  logli1<-cbind(logli1,nlmecpo$BUGSoutput$sims.list$ynewlogli)
}

cpo1<-colMeans(exp(logli1))
lpml1<-sum(log(cpo1))



logli<-NULL

for (i in 1:133) {
  y.3<-y.2[id.2!=i]
  age.3<-age.2[id.2!=i]
  id.4<-id.2[id.2!=i]
  id.3<-id.4
  id.3[id.4>i]<-id.4[id.4>i]-1
  nrec.3<-length(y.3)
  nsub.3<-nsub-1
  ynew<-y.2[id.2==i]
  agenew<-age.2[id.2==i]
  nobs<-length(ynew)
  K=30
  data=list(y=y.3,age=age.3,id=id.3,nrec=nrec.3,
            nsub=nsub.3,ynew=ynew,nobs=nobs,agenew=agenew,K=K)
  
  ##### parameters  
  para<-c("ynewlogli")
  #para<-c("munew","sigmasqnew.i","nu","sigmasq0.i","ynewlogli")

  nlmecpo<-jags(data=data, inits=NULL, parameters.to.save=para, model.file="D:\\dandan\\short course\\examples\\dmd.dir\\nlmeA1_cpo_mean.txt",
                n.chains=1, n.iter=4000, n.burnin=3000,n.thin=1)
  
  logli<-cbind(logli,nlmecpo$BUGSoutput$sims.list$ynewlogli)
}


cpo<-colMeans(exp(logli))
lpml<-sum(log(cpo))

save(logli,file="D:\\dandan\\short course\\examples\\dmd.dir\\nlmecpo.Rdata")
