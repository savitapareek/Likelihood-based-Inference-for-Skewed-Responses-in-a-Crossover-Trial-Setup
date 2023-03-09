## clearing the environment
rm(list=ls())
gc()

## installing the libraries
library(MASS)# for mvrnorm 
library(skewlmm)# for smsn.lmm
library(matlib)# ginv R
library(merDeriv)# for se of varianace components
library(psych) # for trace of matrix
library(fdrtool) # for half normal distribution
library(magic) # for adiag
library(expm)# matrix square root
library(lme4)## for lmer
library(nlme)## for lme fn
library(fBasics)# for colstddev
library(EnvStats)## for qemp
library(mvtnorm)#for dmvnorm
library(matrixcalc)# for vec function
library(heplots)#for cqplot/qqplot
library(sn) ## to generate from skewnormal
library(PDFEstimator)# for estimating the density
library(xtable) ## for latex tables
library(moments) ## for skewness function
library(HLMdiag) # conditional residual
library(ggplot2) # for ggplot
library(reshape2)## for density of each column vector
library(ggpubr)## for ggarrange

## starting time to obtain running time of algorithm
strt=Sys.time()

###*****function for parameter estimation (random error SN)*****
parameter_estimates=function(beth,sigmah,skewparr){
  sigmaeh=sigmah[1]; sigmash=sigmah[2];skewparh=c(skewparr,rep(0,((p*m)-1)))
  delij=skewparh*as.vector((sqrt(1+crossprod(skewparh)))^(-1))
  muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]][[j]]%*%beth))
  dij=sqrt(sigmaeh)*(delij)
  R=(diag(p*m)-crossprod(t(delij)))
  vij=sigmaeh*R+sigmash*crossprod(t(zij))
  tausq=as.numeric(1/(1+(t(dij)%*%solve(vij)%*%dij)))
  eta=lapply(1:s, function(j) sapply(1:no[j], function(i) t(dij)%*%solve(vij)%*%(yij[[j]][i,]-muij[[j]][,i])*tausq))
  T01=lapply(1:s, function(j) sapply(1:no[j], function(i) eta[[j]][i]+dnorm(eta[[j]][i]/sqrt(tausq))/pnorm(eta[[j]][i]/sqrt(tausq))
                                     *sqrt(tausq)))
  T02=lapply(1:s, function(j) sapply(1:no[j], function(i) (eta[[j]][i]**2)+tausq +dnorm(eta[[j]][i]/sqrt(tausq))
                                     /pnorm(eta[[j]][i]/sqrt(tausq))*eta[[j]][i]*sqrt(tausq)))
  
  ## for beta estimates
  bethh=(solve(Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(xij[[j]][[i]])%*%solve(vij)%*%xij[[j]][[i]])))))
         %*%Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(xij[[j]][[i]])%*%solve(vij)%*%(yij[[j]][i,]-dij*T01[[j]][i]))))))
  
  sbet= Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]][[i]])%*%solve(vij)%*%xij[[j]][[i]]))))
  ## parameter estimation of variance components using NR equation
  vo=1
  repeat{
    vo=vo+1
    sigp=c(sigmaeh,sigmash,skewparr);skewparh=c(skewparr,rep(0,(p*m-1)))
    delij=skewparh*as.vector((sqrt(1+crossprod(skewparh)))^(-1))
    muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]][[j]]%*%bethh))
    dij=sqrt(sigmaeh)*(delij)
    R=(diag(p*m)-crossprod(t(delij)))
    vij=sigmaeh*R+sigmash*crossprod(t(zij))
    tausq=as.numeric(1/(1+(t(dij)%*%solve(vij)%*%dij)))
    eta=lapply(1:s, function(j) sapply(1:no[j], function(i) t(dij)%*%solve(vij)%*%(yij[[j]][i,]-muij[[j]][,i])*tausq))
    T01=lapply(1:s, function(j) sapply(1:no[j], function(i) eta[[j]][i]+dnorm(eta[[j]][i]/sqrt(tausq))/pnorm(eta[[j]][i]/sqrt(tausq))
                                       *sqrt(tausq)))
    T02=lapply(1:s, function(j) sapply(1:no[j], function(i) (eta[[j]][i]**2)+tausq +dnorm(eta[[j]][i]/sqrt(tausq))
                                       /pnorm(eta[[j]][i]/sqrt(tausq))*eta[[j]][i]*sqrt(tausq)))
    
    ## derivative of delij and R with respect to \lambda 
    delel=c(1,rep(0,(pm-1)))/(1+skewparr**2)^(3/2); Rl=-2*delij%*%t(delel)
    
    ## derivative of vij with respect to sigma_e^2, sigma_s^2, \lambda
    dvij=list(R,zij%*%t(zij),sigmaeh*Rl)
    
    ## derivative of dij with respect to sigma_e^2, sigma_s^2, \lambda
    ddij=list(delij*(1/(2*sqrt(sigmaeh))),rep(0,p*m),sqrt(sigmaeh)*delel)
    
    ## derivative of solve(vij) with respect to sigma_e^2, sigma_s^2, \lambda
    nopar=length(sigmah)+1## no of parameters to estimate except beta's
    dvijinv=lapply(1:nopar, function(i) -solve(vij)%*%dvij[[i]]%*%solve(vij))
    
    ## derivative of Q function with respect to sigma_e^2, sigma_s^2, \lambda
    dQ_esvrl=sapply(1:nopar, function(k) -.5*(sum(no)*(tr(solve(vij)%*%dvij[[k]]))
                                              +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(yij[[j]][i,]-muij[[j]][,i])%*%dvijinv[[k]]%*%(yij[[j]][i,]-muij[[j]][,i]-2*dij*T01[[j]][i])))))
                                              +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(yij[[j]][i,]-muij[[j]][,i])%*%solve(vij)%*%(-2*ddij[[k]]*T01[[j]][i])))))
                                              +t(dij)%*%dvijinv[[k]]%*%dij*sum(sapply(T02, sum))+2*t(dij)%*%solve(vij)%*%ddij[[k]]*sum(sapply(T02, sum))))
    
    ## second derivative of vij with respect to sigma_e^2, sigma_s^2, \lambda
    delell=-3*skewparr*c(1,rep(0,(pm-1)))*(1+skewparr**2)^(-5/2) ## second derivative of delb
    Rll=-2*delij*delell-2*(delel%*%t(delel))
    vijee=matrix(0,nrow=p*m,ncol=p*m);vijel=Rl
    vijll=sigmaeh*Rll
    ddvij=list(list(vijee,vijee,vijel),
               list(vijee,vijee),
               list(vijll))
    
    ## second derivative of dij with respect to sigma_e^2, sigma_s^2, \lambda
    dijes=rep(0,p*m);dijee=delij*(-1/(4*sigmaeh**(3/2)))
    dijel=delel*(1/(2*sqrt(sigmaeh)))
    dijll=sqrt(sigmaeh)*delell
    dddij=list(list(dijee,dijes,dijel),
               list(dijes,dijes),
               list(dijll))
    
    ## second derivative of solve(vij) with respect to sigma_e^2, sigma_s^2, \lambda
    ddvijinv=  lapply(1:nopar, function(u) lapply(u:nopar, function(u1)
      (-2*solve(vij)%*%dvij[[u]]%*%dvijinv[[u1]]
       -solve(vij)%*%ddvij[[u]][[u1-u+1]]%*%solve(vij))))
    
    ## second derivative of Q function with respect to sigma_e^2, sigma_s^2, \lambda
    ddQ_esvrl=sapply(1:nopar, function(u) sapply(u:nopar, function(u1)
      (-.5*(sum(no)*(tr(solve(vij)%*%ddvij[[u]][[u1-u+1]]+dvijinv[[u1]]%*%dvij[[u]]))
            +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(yij[[j]][i,]-muij[[j]][,i])%*%ddvijinv[[u]][[u1-u+1]]%*%(yij[[j]][i,]-muij[[j]][,i]-2*dij*T01[[j]][i])))))
            +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(yij[[j]][i,]-muij[[j]][,i])%*%dvijinv[[u]]%*%(-2*ddij[[u1]]*T01[[j]][i])))))
            +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(yij[[j]][i,]-muij[[j]][,i])%*%dvijinv[[u1]]%*%(-2*ddij[[u]]*T01[[j]][i])))))
            +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(yij[[j]][i,]-muij[[j]][,i])%*%solve(vij)%*%(-2*dddij[[u]][[u1-u+1]]*T01[[j]][i])))))
            +t(dij)%*%ddvijinv[[u]][[u1-u+1]] %*%dij*sum(sapply(T02, sum))
            +2*t(dij)%*%dvijinv[[u]]%*%ddij[[u1]]*sum(sapply(T02, sum)) 
            +2*t(dij)%*%dvijinv[[u1]]%*%ddij[[u]]*sum(sapply(T02, sum))
            +2*t(ddij[[u1]])%*%solve(vij)%*%ddij[[u]]*sum(sapply(T02, sum))
            +2*t(dij)%*%solve(vij)%*%dddij[[u]][[u1-u+1]] *sum(sapply(T02, sum))  )) ))
    
    H=matrix(nrow=nopar,ncol=nopar)
    for (i in 1:nopar) {
      H[i,i:nopar]=H[i:nopar,i]=ddQ_esvrl[[i]]
    }
    #  if(R(H)<nopar) sigmahh=c(sigp)-Ginv(H)%*%dQ_esvrl
    #if(R(H)==nopar) 
    sigmahh=c(sigp)-solve(H)%*%dQ_esvrl
    print(cbind(sigmahh,c(tsigma,lam)))
    sigmaeh=sigmahh[1];sigmash=sigmahh[2];skewparr=sigmahh[3]
    if(max(abs(sigmahh))>45|vo>200|
       (max(abs(sigmahh-sigp)))<5*10^-3) break 
  }
  
  ## for standard error as diagonal elements of inverse ofInformation matrix
  dQbet_esvrl=sapply(1:nopar, function(k)  -.5*(t(Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i)  t(yij[[j]][i,]-muij[[j]][[i]])%*%(dvijinv[[k]])%*%(-xij[[j]][[i]]))))))
                                                +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]][[i]])%*%dvijinv[[k]]%*%(yij[[j]][i,]-muij[[j]][[i]]-2*dij*T01[[j]][i])) )))
                                                +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]][[i]])%*%solve(vij)%*%(-2*ddij[[k]]*T01[[j]][i])) ))))
  )
  
  I= matrix(nrow=length(test),ncol=length(test))
  pm=(p+t+m-2+1)
  I[1:pm,1:pm]=sbet
  for (i in 1:nopar) {
    I[1:pm,(pm+i)]=I[(pm+i),1:pm]=dQbet_esvrl[,i]
  }
  I[((pm+1):(pm+nopar)),((pm+1):(pm+nopar))]=H
  ## for standard error as diagonal elements of inverse of Information matrix
  # if(R(I)<length(test))  see=sqrt(abs(diag(Ginv(-I))))
  #if(R(I)==length(test))  
  see=sqrt(abs(diag(solve(-I))))
  return(c(bethh,sigmahh,see)) 
  
}

## defining the model specification
p=t=s=3# no of period or time points = no of sequences
m=4# no of genes
no=c(30,30,30) # no. of subjects in each sequence
pm=p*m

### Treatment matrix for X and Z matrix no of rows are p*m : i:no of seq, j= no of subjects in seq i.
treatment=list()
treatment[[1]]=c(1,2,3)
treatment[[2]]=c(2,3,1)
treatment[[3]]=c(3,1,2)
Trt=list()
for (i in 1:s) {
  k=0; k3=m
  Trt[[i]]=matrix(nrow=p*m,ncol = (t-1))
  arr=array(dim = (t-1))
  for (j in 1:p) {
    if(treatment[[i]][j]==1)
      arr[1:(t-1)]=0
    for (tt in 2:t) {
      if(treatment[[i]][j]==tt)
      {arr[tt-1]=1
      arr[-(tt-1)]=0
      }
    }
    Trt[[i]][((k+1):k3),]=kronecker(rep(1,m),rbind(arr))
    k=j*m
    k3=(j+1)*m
  }
}

##individual level co variate #taking values 0, 1 and 2
covv=lapply(1:s, function(i)  c(rep(0,no[i]-s*floor(no[i]/3)),rep(c(0,1,2),each=floor(no[i]/3))))
covv

## X and Z matrix
xij=list()
for (i in 1:s) {
  xij[[i]]=list()
  for (j in 1:no[i]) {
    xij[[i]][[j]]=cbind(rep(1,p*m),rbind(matrix(data=0,nrow = m,ncol = p-1),kronecker(diag(p-1),rep(1,m))),Trt[[i]],kronecker(rep(1,p),rbind(rep(0,m-1),diag(m-1)))
                        ,rep(covv[[i]][j],p*m)  )
  }
}
xij
zij=cbind(rep(1,p*m))

### true values of parameters
cat("enter no of fixed parameters=", (p+t+m-2+1))
#tbet=scan(nmax = (p+t+m-2))
tbet=c(1,2.4, 1.1, .9, 2.1, 1.5, 2, 3.4, 1.8)
sigmas=.8^2; omegsq=2
lam=3
skewpar=c(lam,rep(0,(p*m-1)))
test=c(tbet,omegsq,sigmas,lam);tsigma=c(omegsq,sigmas)

## no of simulations
nr=500

## using model eqn to generate the responses
v1ij=omegsq*diag(p*m)
tdele=skewpar*as.vector(sqrt((1+crossprod(skewpar)))^(-1))
tdij=sqrt(omegsq)*(tdele)
tbij=lapply(1:s, function(i) rnorm(nr*no[i],0,sqrt(sigmas)))
dpM=list(xi=c(0,rep(0,(pm-1))), Omega= v1ij, alpha=skewpar)
teij=lapply(1:s, function(i) rmsn(nr*no[i], dp=dpM))
cpM=dp2cp(dpM,family = "SN")## actual mean, variance and skewness vector of random error
truint=tbet[1]+cpM$mean[1]
test=c(truint,test[-1])
cyij=lapply(1:nr, function(kk) 
  lapply(1:s, function(i) t(sapply(1:no[i], function(j) xij[[i]][[j]]%*%tbet+zij*tbij[[i]][((kk-1)*no[i]+j)]+teij[[i]][((kk-1)*no[i]+j),]))))

###***define matrix to store output***
festi=matrix(nrow=nr,ncol=(2*length(test)))## for proposed algo
festi1=matrix(nrow=nr,ncol=(2*(length(test)-1)))## for lme fitting
aic_bic=matrix(nrow = nr,ncol=4)## proposed and lme fit
## mahalanobis distance and rij vector for each simulated data
maha_dist_n=maha_dist_o=matrix(nrow=nr,ncol=sum(no)) 
resd_n=resd_o=matrix(nrow=nr,ncol=p*m*sum(no))

## simulation code
for (sim in 1:nr) {
  ### variance-covariance matrix of error and random effects
  tteij=lapply(1:s, function(i) teij[[i]][((sim-1)*no[i]+1):(sim*no[i]),])
  yij=cyij[[sim]]
  
  ## data frame for using in lmer function
  resp=unlist(lapply(yij, t))
  period=rep(rep(c("per1","per2","per3"),each=m),sum(no))
  trt=c(rep(rep(c("trt1","trt2","trt3"),each=m),no[1]),
        rep(rep(c("trt2","trt3","trt1"),each=m),no[2]),
        rep(rep(c("trt3","trt1","trt2"),each=m),no[3]))
  gene=rep(rep(paste("gene",1:m),p),sum(no))
  sub=rep(paste("sub",1:sum(no)),each=(p*m))
  cov1=as.vector(sapply(1:s, function(i) rep(covv[[i]],each=p*m)))
  data=data.frame(resp,period,trt,gene,sub,cov1)
  
  ## model fitting based on lme function
  fm1=lme(resp~period+trt+gene+cov1,data=data,random = ~1|sub,method = "ML")
  ## model fitting based lmer function
  model=lmer(resp~period+trt+gene+cov1+ (1|sub),REML = F)# intercept vary randomly with the subjects
  s1=summary(model)
  festi1[sim,]=c(s1$coefficients[,1],as.numeric(VarCorr(fm1)[c(2,1),1]),
                 s1$coefficients[,2],as.numeric(VarCorr(fm1)[c(2,1),2]))
  
  ## aic bic of lme fitting
  tot=p*m*sum(no);nofr=length(test)-1
  aic_bic[sim,1:2]=c( -2*s1$logLik+2*nofr,-2*s1$logLik+log(tot)*nofr )
  
  ## compute mahalanobis distance d_distance
  bett=s1$coefficients[,1]
  sigee=as.numeric(VarCorr(fm1)[c(2),1])
  sigss=as.numeric(VarCorr(fm1)[c(1),1])#;lamm=pp[length(tbet)+3]
  si=sigee*diag(p*m)#;bett=pp[1:length(tbet)]
  sigmai=si+sigss*crossprod(t(zij))
  ddin=unlist(lapply(1:s, function(j) t(sapply(1:no[j], function(i) (t(yij[[j]][i,]-xij[[j]][[i]]%*%bett)%*%solve(sigmai)%*%
                                                                       (yij[[j]][i,]-xij[[j]][[i]]%*%bett))))))
  
  maha_dist_n[sim,]=ddin
  
  ## residuals of lme fitting
  conre=resid_conditional(model)
  resd_n[sim,]=as.numeric(conre)
  
  ###***now using EM algo function with proposed algorithm***&&&%%$$##
  beth=tbet-.1;sigmah=tsigma-.06;skewparr=lam+.05
  ## Now for parameter estimation EM NR method is used.
  pp=parameter_estimates(beth,sigmah,skewparr)
  v0=1
  repeat{
    v0=v0+1
    bethh=pp[1:(p+m+t-2+1)]; sigmahh=pp[(p+t+m-1+1):(p+m+t+1)];skewparh=pp[(p+m+t+1+1)]
    pp=parameter_estimates(bethh,sigmahh,skewparh)
    print(cbind(pp[1:length(test)],test))
    if(max(abs(pp))>40|
       max(abs(pp[1:length(test)]-c(bethh,sigmahh,skewparh)))<5*10**-3) break
  }
  cat("convergence at iteration no",v0)
  la=pp[(length(tbet)+3)]
  delee=la*as.vector((sqrt(1+la^2))^(-1))
  betcp=pp[1]+(sqrt(2/pi)*sqrt(pp[(length(tbet)+1)])*delee)
  festi[sim,]=c(betcp,pp[-1]) ## proposed model estimates
  cbind(test,c(betcp,pp[2:length(test)]),pp[-c(1:length(test))])
  
  # # calculate aic using log likelihood as q fn 
  bett=pp[1:length(tbet)];nofr=length(test)-1
  sigee=pp[length(tbet)+1];sigss=pp[length(tbet)+2];lamm=pp[length(tbet)+3]
  lambr=c(lamm, rep(0,(pm-1)))
  si=sigee*diag(p*m)#;bett=pp[1:length(tbet)]
  sigmai=si+sigss*crossprod(t(zij))
  delb=lambr*as.vector(sqrt((1+crossprod(lambr)))^(-1))
  delbar=sqrt(sigee)*delb;di=delbar
  bigsi=sigmai-crossprod(t(di))# same as vij
  taui=as.numeric(1/(1+t(di)%*%solve(bigsi)%*%di))
  etai=lapply(1:s, function(j) sapply(1:no[j], function(i)
    taui*t(di)%*%solve(bigsi)%*%(yij[[j]][i,]-xij[[j]][[i]]%*%bett)))
  tihat=lapply(1:s, function(j) sapply(1:no[j], function(i)
    etai[[j]][i]+sqrt(taui)*dnorm(etai[[j]][i]/sqrt(taui))/pnorm(etai[[j]][i]/sqrt(taui))))
  tisqhat=lapply(1:s, function(j) sapply(1:no[j], function(i)
    etai[[j]][i]^2+taui+etai[[j]][i]*sqrt(taui)*dnorm(etai[[j]][i]/sqrt(taui))/pnorm(etai[[j]][i]/sqrt(taui))))
  
  qfn= -.5*(sum(no)*(p*m*log(2*pi)+log(det(bigsi))-(log(2)-log(pi)))+Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i)
    t(yij[[j]][i,]-xij[[j]][[i]]%*%bett)%*%solve(bigsi)%*%(yij[[j]][i,]-xij[[j]][[i]]%*%bett-2*tihat[[j]][i]*di)  )))) +
      Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i)
        (tisqhat[[j]][i]+t(di)%*%solve(bigsi)%*%di*tisqhat[[j]][i]))))))
  
  aic_o=-2*qfn+2*nofr
  bic_o=-2*qfn+log(tot)*nofr
  
  aic_bic[sim,3:4]=c(aic_o,bic_o)
  
  ### residual vector based on r_ij
  sn_condresid=lapply(1:s, function(j) (sapply(1:no[j], function(i) (yij[[j]][i,]-xij[[j]][[i]]%*%bett-di*tihat[[j]][i]))))
  resd_o[sim,]=unlist(sn_condresid)
  
  ## compute mahalanobis distance 
  maha_dist_o[sim,]=unlist(lapply(1:s, function(j) t(sapply(1:no[j], function(i) (t(yij[[j]][i,]-xij[[j]][[i]]%*%bett)%*%solve(sigmai)%*%
                                                                                    (yij[[j]][i,]-xij[[j]][[i]]%*%bett))))))
}

## remove the samples which gave bad estimates
sim=sim-1
rm1=which(abs(festi[1:sim,length(test)])>10|(festi[1:sim,length(test)])<.010)
rm2=which(abs(festi[1:sim,(length(test)-1)])>10)
rm3=which((festi[1:sim,(1)])<.10)
rmm=unique(c(rm1,rm2,rm3))
length(rmm)/sim## % of bad samples

## compute the average estimation results based on remaining samples
ns=200## number of simulation you want to have 
aicc=aic_bic[(1:sim)[-rmm][1:ns],]
#lbh_nn=lbh_n[(1:sim)[-rmm][1:ns],];lbh_oo=lbh_o[(1:sim)[-rmm][1:ns],]
festii=festi[(1:sim)[-rmm][1:ns],];festii1=festi1[(1:sim)[-rmm][1:ns],]
maha_dist_nn=colMeans(maha_dist_n[(1:sim)[-rmm][1:ns],])
maha_dist_oo=colMeans(maha_dist_o[(1:sim)[-rmm][1:ns],])
resd_nn=colMeans(resd_n[(1:sim)[-rmm][1:ns],])
resd_oo=colMeans(resd_o[(1:sim)[-rmm][1:ns],])

### qqplot of mahalanobis distance
dataa=data.frame(maha_dist_nn,maha_dist_oo)
colnames(dataa)=c("normal","skew error")
model=rep(c("N fit","SN fit"),each=(sum(no)))
ds1=c(maha_dist_nn,maha_dist_oo)
dff=data.frame(ds1,model)

pl= ggplot(dff, aes(sample = ds1,shape=model,linetype=model))+
  stat_qq(distribution = stats::qchisq,dparams = p*m,size=0.7) + 
  stat_qq_line(distribution = stats::qchisq,dparams = p*m)+
  theme_classic()+
  xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  theme(axis.title=element_text(size=7),plot.title = element_text(size=8))+
  ggtitle("Chi-square Q-Q plot of \n Mahalanobis distances")+
  theme(legend.position = "none")+
  scale_shape_manual(values=c(1,2))
pl

##qq plot of standardized residuals of lme fit and proposed model
std_resd_e=scale(unlist(resd_oo))
std_resd_n=scale(unlist(resd_nn))
df11=data.frame(std_resd_n,std_resd_e)

model=rep(c("N fit","SN fit"),each=(p*m*sum(no)))
ds11=c(std_resd_n,std_resd_e)
dff1=data.frame(ds11,model)

fig= ggplot(dff1, aes(sample = ds11,shape=model,linetype=model))+
  stat_qq(distribution = stats::qnorm,size=0.7) + 
  stat_qq_line(distribution = stats::qnorm)+
  theme_classic()+
  xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  theme(axis.title=element_text(size=7),plot.title = element_text(size=8))+
  ggtitle("Normal Q-Q plot of \n standardized residuals")+
  theme(legend.position = c(.85, .335),legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.key = element_rect(colour = NA, fill = NA))+
  scale_shape_manual(values=c(1,2))
fig 

## arrange both qqplot in one chart
figg=ggarrange(pl,fig)
figg
## final output
fout=cbind(festii,festii1,aicc)
write.csv(fout,"D:\ ranerr out_30.csv")


length(which(aicc[,1]>(aicc[,3]-2*nofr)))/nrow(aicc)*100## AIC glosup formula


## for printing the estimate table
absbias_o= abs(rowMeans(sapply(1:nrow(festii1), function(i) abs(festii[i,1:length(test)]-test))))
absbias_n= abs(rowMeans(sapply(1:nrow(festii1), function(i) abs(festii1[i,1:length(test)]-test))))
snn=cbind(c(test),colMeans(festii[,1:length(test)]),
          colMeans(festii[,-c(1:length(test))]),
          absbias_o)
nn=cbind(colMeans(festii1[,1:(length(test)-1)]),c(test[1:(length(test)-1)]),
         colMeans(festii1[,-c(1:(length(test)-1))]),
         absbias_n[-length(test)])


ll1=cbind(snn,rbind(nn[,-2],0))
ll1
write.csv(ll1,"D:\ errskew2_30.csv")
Sys.time()-strt

#####**** make the empirical densities of ml estimates*****
no_par=12
sim100=festi[(1:sim)[-rmm][1:200],1:no_par]

## first we will work on 100 simulations
esti=sim100

## standardize the estimates
trep=nrow(esti)
obs=matrix(nrow = no_par,ncol=nrow(esti))
for (i in 1:no_par) {
  obs[i,]=(esti[,i]-mean(esti[,i]))/sqrt(var(esti[,i]))
}

## Find the empirical density corresponding to each parameter
x <- qemp(p = seq(0, 1, len = trep), obs[1,]) 
y <- demp(x, obs[1,]) 
x2 <- qemp(p = seq(0, 1, len = trep), obs[2,]) 
y2 <- demp(x2, obs[2,]) 
x3 <- qemp(p = seq(0, 1, len = trep), obs[3,]) 
y3 <- demp(x3, obs[3,]) 
x4 <- qemp(p = seq(0, 1, len = trep), obs[4,]) 
y4 <- demp(x4, obs[4,]) 
x5 <- qemp(p = seq(0, 1, len = trep), obs[5,]) 
y5 <- demp(x5, obs[5,]) 
x6 <- qemp(p = seq(0, 1, len = trep), obs[6,]) 
y6 <- demp(x6, obs[6,]) 
x7 <- qemp(p = seq(0, 1, len = trep), obs[7,]) 
y7 <- demp(x7, obs[7,]) 
x8 <- qemp(p = seq(0, 1, len = trep), obs[8,]) 
y8 <- demp(x8, obs[8,]) 
x9 <- qemp(p = seq(0, 1, len = trep), obs[9,]) 
y9 <- demp(x8, obs[9,]) 
x10 <- qemp(p = seq(0, 1, len = trep), obs[10,]) 
y10<- demp(x8, obs[10,]) 
x11 <- qemp(p = seq(0, 1, len = trep), obs[11,]) 
y11<- demp(x8, obs[11,]) 
x12 <- qemp(p = seq(0, 1, len = trep), obs[12,]) 
y12<- demp(x8, obs[12,]) 
## data frame and ggplot graph
df33 <-data.frame(parameter=rep(c("Intercept","Period2","Period3",
                                  "Trt2","Trt3","Gene2","Gene3","Gene4","ind_level_covariate"), each=trep),
                  Parameter_value=c(x,x2,x3,x4,x5,x6,x7,x8,x9),
                  Relative_frequency=c(y,y2,y3,y4,y5,y6,y7,y8,y9))
p1=ggplot(df33, aes(Parameter_value,Relative_frequency,group=parameter)) +
  geom_line(aes(linetype=parameter),lwd=.09)+ theme_classic()+ 
  ggtitle("Empirical density plot of \n maximum likelihood estimates")+
  theme(axis.title=element_text(size=7),plot.title = element_text(size=8))+
  theme(legend.position = "none")
p1

ggarrange(figg, p1,nrow=2)

