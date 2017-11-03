
library(R2jags)
library(coda)
library(lattice)



# 1. Read observed data
### Change the directory
GH<-read.table("D://dandan//short course//examples//GH//GH.txt", header=T, sep=' ')  # data frame containing the observed GH data
N<-nrow(GH)  							# N: total number of observations
yall<-as.matrix(GH[,1:3])						# N by 3 matrix each comlumn of which corresponds to one
# of the three observation times
trtall<-as.integer(GH[,4])						# numeric vector indicating the assigned treatment to each subject

# RUN BNP MODEL 1 -- treatment 1
# 2. Reorder the data for treatment 1
yall2<-yall[trtall==1,]
# standadize the data
yall1<-(yall2-mean(c(yall2),na.rm=TRUE))/sd(c(yall2),na.rm=TRUE)*sqrt(0.5)
ymean<-mean(c(yall2),na.rm=TRUE)
ysd<-sd(c(yall2),na.rm=TRUE)/sqrt(0.5)


y1<-yall1[is.na(yall1[,2])==TRUE&is.na(yall1[,3])==TRUE,]
y2<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==TRUE,]
y3<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==FALSE,]
y<-rbind(y1,y2,y3)
N1<-nrow(y1)
N2<-nrow(rbind(y1,y2))
N<-nrow(y)
K<-50
g1<-(sd(yall2[is.na(yall2[,1])==FALSE,1])/ysd)^2
g2<-(summary(lm(yall2[,2]~yall2[,1]))$sigma/ysd)^2
g3<-(summary(lm(yall2[,3]~yall2[,1]+yall2[,2]))$sigma/ysd)^2

# 3. run the model for treatment 1
param<-c('cluster','mmu1','mmu2','mmu3','logli')  		
data=list(y=y,N=N,K=K,N1=N1,N2=N2,g1=g1,g2=g2,g3=g3,ymean=ymean,ysd=ysd)
fits1<-jags(data=data,  parameters.to.save=param, 
            model.file="D://dandan//short course//examples//GH//DPM//DPMmod12.txt",
            n.chains=1, n.iter=15000, n.burnin=5000,n.thin=1)


# RUN BNP MODEL 1 -- treatment 2
# 4. Reorder the data for treatment 2
yall2<-yall[trtall==2,]
# standadize the data
yall1<-(yall2-mean(c(yall2),na.rm=TRUE))/sd(c(yall2),na.rm=TRUE)*sqrt(0.5)
ymean<-mean(c(yall2),na.rm=TRUE)
ysd<-sd(c(yall2),na.rm=TRUE)/sqrt(0.5)

y1<-yall1[is.na(yall1[,2])==TRUE&is.na(yall1[,3])==TRUE,]
y2<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==TRUE,]
y3<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==FALSE,]
y<-rbind(y1,y2,y3)
N1<-nrow(y1)
N2<-nrow(rbind(y1,y2))
N<-nrow(y)
K<-50
g1<-(sd(yall2[is.na(yall2[,1])==FALSE,1])/ysd)^2
g2<-(summary(lm(yall2[,2]~yall2[,1]))$sigma/ysd)^2
g3<-(summary(lm(yall2[,3]~yall2[,1]+yall2[,2]))$sigma/ysd)^2

# 5. run the model for treatment 2
data=list(y=y,N=N,K=K,N1=N1,N2=N2,g1=g1,g2=g2,g3=g3,ymean=ymean,ysd=ysd)
fits2<-jags(data=data,  parameters.to.save=param, 
            model.file="D://dandan//short course//examples//GH//DPM//DPMmod12.txt",
            n.chains=1, n.iter=15000, n.burnin=5000,n.thin=1)


# RUN BNP MODEL 1 -- treatment 3
# 6. Reorder the data for treatment 3
yall2<-yall[trtall==3,]
# standadize the data
yall1<-(yall2-mean(c(yall2),na.rm=TRUE))/sd(c(yall2),na.rm=TRUE)*sqrt(0.5)
ymean<-mean(c(yall2),na.rm=TRUE)
ysd<-sd(c(yall2),na.rm=TRUE)/sqrt(0.5)

y1<-yall1[is.na(yall1[,2])==TRUE&is.na(yall1[,3])==TRUE,]
y2<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==TRUE,]
y3<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==FALSE,]
y<-rbind(y1,y2,y3)
N1<-nrow(y1)
N2<-nrow(rbind(y1,y2))
N<-nrow(y)
K<-50
g1<-(sd(yall2[is.na(yall2[,1])==FALSE,1])/ysd)^2
g2<-(summary(lm(yall2[,2]~yall2[,1]))$sigma/ysd)^2
g3<-(summary(lm(yall2[,3]~yall2[,1]+yall2[,2]))$sigma/ysd)^2


# 7. run the model for treatment 3
data=list(y=y,N=N,K=K,N1=N1,N2=N2,g1=g1,g2=g2,g3=g3,ymean=ymean,ysd=ysd)
fits3<-jags(data=data,  parameters.to.save=param, 
            model.file="D://dandan//short course//examples//GH//DPM//DPMmod12.txt",
            n.chains=1, n.iter=15000, n.burnin=5000,n.thin=1)


# RUN BNP MODEL 1 -- treatment 4
# 8. Reorder the data for treatment 4
yall2<-yall[trtall==4,]
# standadize the data
yall1<-(yall2-mean(c(yall2),na.rm=TRUE))/sd(c(yall2),na.rm=TRUE)*sqrt(0.5)
ymean<-mean(c(yall2),na.rm=TRUE)
ysd<-sd(c(yall2),na.rm=TRUE)/sqrt(0.5)

y1<-yall1[is.na(yall1[,2])==TRUE&is.na(yall1[,3])==TRUE,]
y2<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==TRUE,]
y3<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==FALSE,]
y<-rbind(y1,y2,y3)
N1<-nrow(y1)
N2<-nrow(rbind(y1,y2))
N<-nrow(y)
K<-50
g1<-(sd(yall2[is.na(yall2[,1])==FALSE,1])/ysd)^2
g2<-(summary(lm(yall2[,2]~yall2[,1]))$sigma/ysd)^2
g3<-(summary(lm(yall2[,3]~yall2[,1]+yall2[,2]))$sigma/ysd)^2


# 9. run the model for treatment 4
data=list(y=y,N=N,K=K,N1=N1,N2=N2,g1=g1,g2=g2,g3=g3,ymean=ymean,ysd=ysd)
fits4<-jags(data=data,  parameters.to.save=param, 
            model.file="D://dandan//short course//examples//GH//DPM//DPMmod12.txt",
            n.chains=1, n.iter=15000, n.burnin=5000,n.thin=1)

LPML4<-sum(-log(colMeans(1/exp(fits1$BUGSoutput$sims.list$logli))))+
  sum(-log(colMeans(1/exp(fits2$BUGSoutput$sims.list$logli))))+
  sum(-log(colMeans(1/exp(fits3$BUGSoutput$sims.list$logli))))+
  sum(-log(colMeans(1/exp(fits4$BUGSoutput$sims.list$logli))))


result1<-fits1$BUGSoutput$sims.matrix[,41:43]
colMeans(result1)
apply(result1,2,function (x) {quantile(x,0.025)})
apply(result1,2,function (x) {quantile(x,0.975)})

result2<-fits2$BUGSoutput$sims.matrix[,44:46]
colMeans(result2)
apply(result2,2,function (x) {quantile(x,0.025)})
apply(result2,2,function (x) {quantile(x,0.975)})

result3<-fits3$BUGSoutput$sims.matrix[,43:45]
colMeans(result3)
apply(result3,2,function (x) {quantile(x,0.025)})
apply(result3,2,function (x) {quantile(x,0.975)})

result4<-fits4$BUGSoutput$sims.matrix[,44:46]
colMeans(result4)
apply(result4,2,function (x) {quantile(x,0.025)})
apply(result4,2,function (x) {quantile(x,0.975)})



#  run the BNP model 2
yall2<-yall
yall1<-(yall2-mean(c(yall2),na.rm=TRUE))/sd(c(yall2),na.rm=TRUE)*sqrt(0.5)
ymean<-mean(c(yall2),na.rm=TRUE)
ysd<-sd(c(yall2),na.rm=TRUE)/sqrt(0.5)

y1<-yall1[is.na(yall1[,2])==TRUE&is.na(yall1[,3])==TRUE,]
y2<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==TRUE,]
y3<-yall1[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==FALSE,]
y<-rbind(y1,y2,y3)
N1<-nrow(y1)
N2<-nrow(rbind(y1,y2))
N<-nrow(y)
trt<-c(trtall[is.na(yall1[,2])==TRUE&is.na(yall1[,3])==TRUE],
       trtall[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==TRUE],
       trtall[is.na(yall1[,2])==FALSE&is.na(yall1[,3])==FALSE])
K<-50
g1<-(sd(yall2[is.na(yall2[,1])==FALSE,1])/ysd)^2
g2<-(summary(lm(yall2[,2]~yall2[,1]))$sigma/ysd)^2
g3<-(summary(lm(yall2[,3]~yall2[,1]+yall2[,2]))$sigma/ysd)^2


param<-c('cluster','mmu1','mmu2','mmu3','logli')    	
data=list(y=y,N=N,K=K,N1=N1,N2=N2,g1=g1,g2=g2,g3=g3,ymean=ymean,ysd=ysd,trt=trt)
fits5<-jags(data=data,  parameters.to.save=param, 
            model.file="D://dandan//short course//examples//GH//DPM//DPMmod16.txt",
            n.chains=1, n.iter=15000, n.burnin=5000,n.thin=1)


LPML5<-sum(-log(colMeans(1/exp(fits5$BUGSoutput$sims.list$logli))))

result5<-fits5$BUGSoutput$sims.matrix[,163:174]
colMeans(result5)
apply(result5,2,function (x) {quantile(x,0.025)})
apply(result5,2,function (x) {quantile(x,0.975)})


# run BNP model 3
param<-c('cluster','mmu1','mmu2','mmu3','logli')    	
data=list(y=y,N=N,K=K,N1=N1,N2=N2,g1=g1,g2=g2,g3=g3,ymean=ymean,ysd=ysd,trt=trt)
fits6<-jags(data=data,  parameters.to.save=param, 
            model.file="D://dandan//short course//examples//GH//DPM//DPMmod17.txt",
            n.chains=1, n.iter=15000, n.burnin=5000,n.thin=1)


LPML6<-sum(-log(colMeans(1/exp(fits6$BUGSoutput$sims.list$logli))))

result6<-fits6$BUGSoutput$sims.matrix[,163:174]
colMeans(result6)
apply(result6,2,function (x) {quantile(x,0.025)})
apply(result6,2,function (x) {quantile(x,0.975)})


table(trtall)







