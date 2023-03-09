## clean the environment
rm(list=ls())

## install the relevant libraries
library(moments)
library(lme4)
library(MASS)## for boxcox
library(nlme)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rcompanion)
library(dgof)# for ks.test
library(tidyverse)
library(Metrics)# for relative error
library(moments)## for skewness
library(skewlmm)# for smsn.lmm fn
library(HLMdiag)# for marginal residuals and conditional resiuals

## read the data
dat=read.csv(file.choose())## read anova 3 way all seq excel from csv files from vol D
head(dat)
subb=paste("Subject",1:12)
Subject=factor(rep(paste("S",1:12, sep=""),each=30))
Sequence=factor(rep(c("Sequence 1","Sequence 2","Sequence 3"),each=120))
Treatment=factor(c(rep(rep(c("T2","T1","T3"),each=10),4),
       rep(rep(c("T3","T2","T1"),each=10),4),
       rep(rep(c("T1","T3","T2"),each=10),4)))
Gene=factor(rep(paste("Gene",1:10,sep=""),36))
Period=factor(rep(rep(paste("P",1:3,sep=""),each=10),12))
Response=dat$Response_allseq

## data frame having all the variables period, treatment gene and responses
data=data.frame(Response,Subject,Sequence,Treatment,Gene,Period)
attach(data)
head(data)
View(data)
str(data)

### histogram and qqplot of original response
df=data
p1=ggplot(df, aes(Response)) +        # Draw density
  geom_density() +
  stat_function(fun = dnorm,
                args = list(mean = mean(Response),
                            sd = sd(Response)),
                col = "black",linetype = "dashed",
                size = .5)+theme_classic()
p11=p1+ylab("Density")+xlab("Response \n mRNA expression levels (pg/ml)")+
  theme(axis.title=element_text(size=10))
p11
## draw qqplot
p <- ggplot(df, aes(sample = Response))+ stat_qq(shape=1) + stat_qq_line()+theme_classic()
p2=p+xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  theme(axis.title=element_text(size=10))

figure <- ggarrange(p11, p2,
                    
                    ncol = 2, nrow = 1,labels = c("A","B"),
                    font.label=list(color="black",size=9))
figure

## test using ks.yest and shapiro wilk test
shapiro.test(Response) # p-value<0.0001 so data is not normal
ks.test(Response,"pnorm",mean=mean(Response), sd=sd(Response))

## now do the box-cox transformation
rr=boxcox(Response~Gene+Period+Treatment)
rr1=boxcox(Response~1)
rr1$x[which.max(rr1$y)]
vaa=rr$x[which.max(rr$y)]
qqnorm(Response^vaa)
shapiro.test(Response^vaa) ## box-cox transformation also suggest non-normality

##or sqrt transformation
shapiro.test(sqrt(Response)) ## sqrt transformation is also non-normal

## test using ks.yest and shapiro wilk test for log response
shapiro.test(log(Response))
ks.test(log(Response),"pnorm",mean=mean(log(Response)), sd=sd(log(Response)))
skewness((Response))
  
## desnity and qqplot of log responses
p1=ggplot(df, aes(log(Response))) +        # Draw density
  geom_density() +
  stat_function(fun = dnorm,
                args = list(mean = mean(log(Response)),
                            sd = sd(log(Response))),
                col = "black",linetype="dashed",
                size = .5)+theme_classic()
p11=p1+ylab("Density")+xlab("log(Response) \n log(mRNA expression levels (pg/ml))")+
  theme(axis.title=element_text(size=10))
p <- ggplot(df, aes(sample = log(Response)))+ stat_qq(shape=1) + stat_qq_line()+theme_classic()
p2=p+xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  theme(axis.title=element_text(size=10))

figure1 <- ggarrange(p11, p2,
                    
                    ncol = 2, nrow = 1,labels = c("C","D"),
                    font.label=list(color="black",size=9))
figure1
ggarrange(figure,figure1,nrow=2,ncol=1)


### now plot the effect of period and response for each subject or subject by period interaction
## period by gene
pergen=ggplot(df, aes(Period,Response,group=Gene,color=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  stat_summary(fun = mean, geom = "point",size=.41,color="black")+
  font("y.text", size = 6)+font("x.text", size = 6)+
  ylab("Response (mRNA expression levels (pg/ml))")+
  theme(axis.title=element_text(size=8))+
  scale_x_discrete(labels = c("P1" = expression("P"[1]),
                              "P2" = expression("P"[2]),
                              "P3" = expression("P"[3])))

pergen
## treatment by gene
trtgen=ggplot(df, aes(Treatment,Response,group=Gene,color=Gene)) +
  stat_summary(fun = mean, geom = "line")+
  stat_summary(fun = mean, geom = "point",size=.41,color="black")+
  font("y.text", size = 6)+font("x.text", size = 6)+
  ylab("Response (mRNA expression levels (pg/ml))")+
  theme(axis.title=element_text(size=8))+
  scale_x_discrete(labels = c("T1" = expression("T"[1]),
                              "T2" = expression("T"[2]),
                              "T3" = expression("T"[3])))

trtgen
##Subject by gene
Subjects=factor(Subject,labels = paste("S",1:12,sep=""))
ndf=data.frame(df,Subjects)
sdf=filter(ndf,Gene==c("Gene10","Gene2","Gene9"))
subgen=ggplot(ndf, aes(Subjects,Response,group=Gene,color=Gene)) +
  stat_summary(fun = mean, geom = "line")+ 
  stat_summary(fun = mean, geom = "point",size=.31,color="black")+
  font("y.text", size = 6)+font("x.text", size = 6)+
  ylab("Response (mRNA expression levels (pg/ml))")+
  theme(axis.title=element_text(size=8),legend.text=element_text(size=8))+
  scale_x_discrete(labels = c("S1" = expression("S"[1]),"S2" = expression("S"[2]),"S3" = expression("S"[3]),
                              "S4" = expression("S"[4]),"S5" = expression("S"[5]),"S6" = expression("S"[6]),
                              "S7" = expression("S"[7]),"S8" = expression("S"[8]),"S9" = expression("S"[9]),
                              "S10" = expression("S"[10]),"S11" = expression("S"[11]),"S12" = expression("S"[12])
                              ))
subgen

## final interaction plot
ggarrange(pergen,trtgen+rremove("ylab"),subgen+rremove("ylab"),nrow = 1,ncol = 3,legend = "bottom",common.legend = T
          ,labels = c("A","B","C"),
          font.label=list(color="black",size=9))#

### fitting the lme, lmer  and smsn.lmm model for original responses&&&&*****%%%$$$
model=lmer((Response)~ Period+ Treatment+Gene+ (1|Subject),data=df,REML = F)# intercept vary randomly with the subjects
model
fm2=lme((Response)~ Period+ Treatment+Gene,data=df,random = ~1|Subject,method = "ML")
fm2
plot(fm2)

fm1 <- smsn.lmm(df,formFixed=(Response)~ Period+ Treatment+Gene,
                groupVar="Subject")
summary(fm1)
plot(fm1)
healy.plot(fm1)

## check the estimated random effects
ar=random.effects(fm2)
qq= as.numeric(substr(row.names(ar),2,3))
rt=sapply(1:length(unique(Subject)), function(i) which(qq==i))
esran= ar[[1]][rt]
dff=data.frame(esran,c(1:12))

p1=ggplot(dff, aes(esran)) +        # Draw density of estimated random intercepts
  geom_density() +
  stat_function(fun = dnorm,
                args = list(mean = mean(dff$esran),
                            sd = sd(dff$esran)),
                col = "black",linetype="dashed",
                size = .5)+theme_classic()+xlab("Empirical estimates of intercepts")
p12=p1+ylab("Density")+xlim(-.1,.1)+
  theme(axis.title=element_text(size=10))
p12
fig1 <- ggplot(dff, aes(sample = esran))+ stat_qq(shape=1) + stat_qq_line()+theme_classic()+
  xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  theme(axis.title=element_text(size=10))
fig1
figure1 <- ggarrange(p12, fig1,

                     ncol = 2, nrow = 1,labels = c("A","B"),
                     font.label=list(color="black",size=9))
figure1
## shapiro wilk test for random effects
shapiro.test(esran)

## now for residuals normality
conre=resid_conditional(model,"studentized")
qqnorm(conre)
hist(conre)
residuals(model)
rstudent(model)
cbind(residuals(model),conre)
boxplot(scale(residuals(fm2)))$out
#scale(residuals(fm2))[262]## 262 is outlying observation so we remove it
df1=data.frame(fitted.values(fm2)[-262],scale(residuals(fm2))[-262])#conre)#
colnames(df1)=c("Fitted values","Standardized residuals")

attach(df1)
fig2= ggplot(df1, aes(x = `Fitted values`, y =`Standardized residuals`)) + geom_point(shape=1)+
geom_hline(yintercept = 0)+theme_classic()+
  theme(axis.title=element_text(size=10))
fig2

fig4= ggplot(df1, aes(`Standardized residuals`)) + geom_density()+
               stat_function(fun = dnorm,
                             args = list(mean = mean(`Standardized residuals`),
                                         sd = sd(`Standardized residuals`)),
                             col = "black",linetype="dashed",
                             size = .5)+theme_classic()
fig4

fig3=ggplot(df1, aes(sample = `Standardized residuals`)) + stat_qq(shape=1) +stat_qq_line()+
  theme_classic()+
  xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  theme(axis.title=element_text(size=10))
fig3

ff2=ggarrange(fig2, fig3,ncol = 2, nrow = 1,labels = c("C","D"),
              font.label=list(color="black",size=9))
ggarrange(figure1,ff2,nrow = 2)

# fpp=ggarrange(fig4,fig3,ncol=2)
# ggarrange(fpp,fig2,nrow=2)

## shapiro wilk test for residuals
shapiro.test(scale(residuals(fm2)))
shapiro.test(residuals(model))

##$$$******%%%%%%%%%%%%%%% to make the healy plot$$$$%%%%@@@####!!!!
healy.plot(fm1)
s=3#no of sequences
p=t=3#no of periods and no of treatments
m=10#no of  genes
no=array()
pm=p*m
no=c(4,4,4)
treatment1=list()
treatment1[[1]]=c(2,1,3)#,3)
treatment1[[2]]=c(3,2,1)#,1)
treatment1[[3]]=c(1,3,2)

######################### now u can run your code at once
Trt=list()
for (i in 1:s) {
  k=0; k3=m
  Trt[[i]]=matrix(nrow=p*m,ncol = (t-1))
  arr=array(dim = t-1)
  for (j in 1:p) {
    if(treatment1[[i]][j]==1)
      arr[1:(t-1)]=0
    for (tt in 2:t) {
      if(treatment1[[i]][j]==tt)
      {arr[tt-1]=1
      arr[-(tt-1)]=0
      }
    }
    Trt[[i]][((k+1):k3),]=kronecker(rep(1,m),rbind(arr))
    k=j*m
    k3=(j+1)*m
  }
}


##############
x1k=list()
for (i in 1:s) {
  x1k[[i]]=cbind(rep(1,p*m),rbind(matrix(data=0,nrow = m,ncol = p-1),kronecker(diag(p-1),rep(1,m))),Trt[[i]],kronecker(rep(1,p),rbind(rep(0,m-1),diag(m-1))))#,carry[[i]]
}

#making z matrix
z1k=cbind(kronecker(rep(1,p),rep(1,m)))
xij=x1k;zij=z1k
y=list()
y[[1]]=Response[(1:120)]
y[[2]]=Response[(121:240)]
y[[3]]=Response[(241:360)]
yij=lapply(1:s, function(i) matrix(y[[i]],ncol=30,byrow=T))
skewparre=fm1$estimates$lambda
bethhe=fm1$estimates$beta[c(1:5,7:14,6)]
sigmashe=as.numeric(fm1$estimates$D)
sigmaehe=fm1$estimates$sigma2
delb=skewparre*as.vector((sqrt(1+skewparre^2))^(-1))
muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]]%*%bethhe))
dij=sqrt(sigmashe)*(delb)*zij
R=(1-delb^2)
vij=sigmaehe*diag(p*m)+R*sigmashe*crossprod(t(zij))
dist=unlist(lapply(1:s, function(j) sapply(1:no[j], function(i) (t(yij[[j]][i,]-xij[[j]]%*%bethhe)%*%vij%*%
                                                                   (yij[[j]][i,]-xij[[j]]%*%bethhe)))))


orobsmahadist_y=sort(qchisq(dist,30))
xval=seq(from=(1/12),to=(12/12),length.out=12)
plot(xval,orobsmahadist_y)
abline(0,0)