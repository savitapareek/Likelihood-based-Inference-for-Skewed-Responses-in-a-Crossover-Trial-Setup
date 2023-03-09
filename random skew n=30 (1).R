## clearing the enviornment
rm(list=ls())
gc()

## installing the libraries
library(MASS)# for mvrnorm 
library(heplots)# for cqplots fun
library(psych) # for trace of matrix
library(fdrtool) # for half normal distribution
library(magic) # for adiag
library(expm)# matrix square root
library(matrixcalc)# for vec function
library(sn) ## to generate from skewnormal
library(xtable) ## for latex tables
library(moments)## for skewness
library(nlme)## for lme fn
library(lme4)## for lmer
library(matlib)# ginv R
library(EnvStats)## for qemp fn
library(lmeInfo)# for lme goodness of fit
library(ggplot2)## for ggplots
library(ggpubr)## for ggarrange
library(resample)#colstdevs
library(skewlmm)# for smsn.lmm
library(PDFEstimator)## for estimating pdf

## starting time of the program
strt=Sys.time()

## parameter estimate using proposed algorithm
parameter_estimates=function(beth,sigmah,skewparr){
  sigmaeh=sigmah[1]; sigmash=sigmah[2]
  delb=skewparr*as.vector((sqrt(1+skewparr^2))^(-1))
  muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]][[j]]%*%beth))
  dij=sqrt(sigmash)*(delb)*zij
  R=(1-delb^2)
  vij=sigmaeh*diag(p*m)+R*sigmash*crossprod(t(zij))
  tausq=as.numeric(1/(1+(t(dij)%*%solve(vij)%*%dij)))
  eta=lapply(1:s, function(j) sapply(1:no[j], function(i) t(dij)%*%solve(vij)%*%(yij[[j]][i,]-muij[[j]][,i])*tausq))
  T01=lapply(1:s, function(j) sapply(1:no[j], function(i) eta[[j]][i]+dnorm(eta[[j]][i]/sqrt(tausq))/pnorm(eta[[j]][i]/sqrt(tausq))
                                     *sqrt(tausq)))
  T02=lapply(1:s, function(j) sapply(1:no[j], function(i) (eta[[j]][i]**2)+tausq +dnorm(eta[[j]][i]/sqrt(tausq))
                                     /pnorm(eta[[j]][i]/sqrt(tausq))*eta[[j]][i]*sqrt(tausq)))
  
  ## for beta estimates
  bethh=(solve(Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(xij[[j]][[i]])%*%solve(vij)%*%xij[[j]][[i]])))))
         %*%Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(xij[[j]][[i]])%*%solve(vij)%*%(yij[[j]][i,]-dij*T01[[j]][i]))))))
  
  ## for information matrix of beta
  sbet= Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]][[i]])%*%solve(vij)%*%xij[[j]][[i]]))))
  
  ## use NR method for estimating variance and skewness parameters
  vo=1
  repeat{
    vo=vo+1
    sigp=c(sigmaeh,sigmash,skewparr)
    delb=skewparr*as.vector((sqrt(1+skewparr^2))^(-1))
    muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]][[j]]%*%bethh))
    dij=sqrt(sigmash)*(delb)*zij
    R=(1-delb^2)
    vij=sigmaeh*diag(p*m)+R*sigmash*crossprod(t(zij))
    tausq=as.numeric(1/(1+(t(dij)%*%solve(vij)%*%dij)))
    eta=lapply(1:s, function(j) sapply(1:no[j], function(i) t(dij)%*%solve(vij)%*%(yij[[j]][i,]-muij[[j]][,i])*tausq))
    T01=lapply(1:s, function(j) sapply(1:no[j], function(i) eta[[j]][i]+dnorm(eta[[j]][i]/sqrt(tausq))/pnorm(eta[[j]][i]/sqrt(tausq))
                                       *sqrt(tausq)))
    T02=lapply(1:s, function(j) sapply(1:no[j], function(i) (eta[[j]][i]**2)+tausq +dnorm(eta[[j]][i]/sqrt(tausq))
                                       /pnorm(eta[[j]][i]/sqrt(tausq))*eta[[j]][i]*sqrt(tausq)))
    
    delbl=1/((1+skewparr**2)^(3/2)); Rl=-2*delb*delbl
    ## derivative of vij with respect to sigma_e^2, sigma_s^2, \lambda
    dvij=list(diag(p*m),R*zij%*%t(zij),sigmash*Rl*zij%*%t(zij))
    
    ## derivative of dij with respect to sigma_e^2, sigma_s^2, \lambda
    ddij=list(rep(0,p*m),zij*delb*(1/(2*sqrt(sigmash))),zij*sqrt(sigmash)*delbl)
    
    ## derivative of solve(vij) with respect to sigma_e^2, sigma_s^2, \lambda
    nopar=length(sigmah)+1## no of parameters to estimate except beta's
    dvijinv=lapply(1:nopar, function(i) -solve(vij)%*%dvij[[i]]%*%solve(vij))
    
    ## derivative of Q function with respect to sigma_e^2, sigma_s^2,  \lambda
    dQ_esvrl=sapply(1:nopar, function(k) -.5*(sum(no)*(tr(solve(vij)%*%dvij[[k]]))
                                              +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(yij[[j]][i,]-muij[[j]][,i])%*%dvijinv[[k]]%*%(yij[[j]][i,]-muij[[j]][,i]-2*dij*T01[[j]][i])))))
                                              +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(yij[[j]][i,]-muij[[j]][,i])%*%solve(vij)%*%(-2*ddij[[k]]*T01[[j]][i])))))
                                              +t(dij)%*%dvijinv[[k]]%*%dij*sum(sapply(T02, sum))+2*t(dij)%*%solve(vij)%*%ddij[[k]]*sum(sapply(T02, sum))))
    
    ## second derivative of vij with respect to sigma_e^2, sigma_s^2, \lambda
    delbll=-3*skewparr*(1+skewparr**2)^(-5/2) ## second derivative of delb
    Rll=-2*delb*delbll-2*(delbl**2)
    vijee=matrix(0,nrow=p*m,ncol=p*m);vijsl=Rl*zij%*%t(zij)
    vijll=sigmash*Rll*zij%*%t(zij)
    ddvij=list(list(vijee,vijee,vijee),
               list(vijee,vijsl),
               list(vijll))
    
    ## second derivative of dij with respect to sigma_e^2, sigma_s^2, \lambda
    dijee=rep(0,p*m);dijss=zij*delb*(-1/(4*sigmash**(3/2)))
    dijsl=zij*delbl*(1/(2*sqrt(sigmash)))
    dijll=zij*sqrt(sigmash)*delbll
    dddij=list(list(dijee,dijee,dijee),
               list(dijss,dijsl),
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
    #if(R(H)<nopar) sigmahh=c(sigp)-Ginv(H)%*%dQ_esvrl
    #if(R(H)==nopar) 
    sigmahh=c(sigp)-solve(H)%*%dQ_esvrl
    sigmaeh=sigmahh[1];sigmash=sigmahh[2];skewparr=sigmahh[3]
    if(max(abs(sigmahh))>50|vo>200|(max(abs(sigmahh-sigp)))<5*10^-3) break
  }
  cat("convergence of variance component, NR eqn", vo)
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
  #if(R(I)<length(test))  see=sqrt(abs(diag(Ginv(-I))))
  #if(R(I)==length(test))  
  see=sqrt(abs(diag(solve(-I))))
  return(c(bethh,sigmahh,see))
}

## defining the model specifications
p=t=s=3# no of period or time points = no of sequences
m=4# no of genes
no=c(30,30,30) # no. of subjects in each sequence
pm=p*m

### X and Z matrix no of rows are p*m : i:no of seq, j= no of subjects in seq i.
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

##individual level co variate
covv=lapply(1:s, function(i)  c(rep(0,no[i]-s*floor(no[i]/3)),rep(c(0,1,2),each=floor(no[i]/3))))
covv

## x and z matrix
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
tbet=c(2,2.4, 1.1, .9, 2.1, 1.5, 2, 3.4, 1.8)
sigmaesq=.85^2; omegsq=3
lam=4
test=c(tbet,sigmaesq,omegsq,lam);tsigma=c(sigmaesq,omegsq)

nr=500 # number of simulations

## using sn function to generate the random effects skew normal variate
tdelb=lam*as.vector((sqrt(1+lam^2))^(-1))
tdij=sqrt(omegsq)*(tdelb)*zij
tv1ij=sigmaesq*diag(p*m)
eij=lapply(1:s, function(i) mvrnorm(nr*no[i],rep(0,pm),tv1ij))
dpM <- c(0,sqrt(omegsq),lam)
cpM=dp2cp(dpM,family = "SN")## actual mean and variances of SN 
truint=tbet[1]+cpM[1]
test=c(truint,test[-1])
tbij= lapply(1:s, function(i) rsn(nr*no[i],dp=dpM))

#and matrix to store output
festi=festi2=matrix(nrow=nr,ncol=(2*length(test)))## for proposed and smsn fit
festi1=matrix(nrow=nr,ncol=(2*(length(test)-1))) ## for lme fit
aic_bic=matrix(nrow = nr,ncol=6)## proposed, smsn and lme fit
## mahalanobis distance and rij vector for each simulated data
maha_dist_n=maha_dist_o=maha_dist_smsn= matrix(nrow=nr,ncol=sum(no)) 
resd_n=resd_o=resd_smsn=matrix(nrow=nr,ncol=p*m*sum(no))

## code for each simulation run
for (sim in 1:nr) {
  ## data generation, read one sample from the above generated data
  tbbij=lapply(1:s, function(i) tbij[[i]][((sim-1)*no[i]+1):(sim*no[i])])
  teij=lapply(1:s, function(i) eij[[i]][((sim-1)*no[i]+1):(sim*no[i]),])
  yij= lapply(1:s, function(i) t(sapply(1:no[i], function(j) xij[[i]][[j]]%*%tbet+zij*tbbij[[i]][j]+teij[[i]][j,])))
  
  ## create the data frame first for lme and smsn.lmm package
  period=rep(rep(c("p1","p2","p3"),each=m),sum(no))
  tyij=lapply(yij, t)
  resp=unlist(tyij)#as.vector(t(yij[[1]]))#
  trt=c(rep(rep(c("t1","t2","t3"),each=m),(no[1])),
        rep(rep(c("t2","t3","t1"),each=m),(no[2])),
        rep(rep(c("t3","t1","t2"),each=m),(no[3])))
  gene=rep(rep(paste("gene",1:m),p),sum(no))
  sub=rep(paste("sub",1:sum(no)),each=(p*m))
  cov1=as.vector(sapply(1:s, function(i) rep(covv[[i]],each=p*m)))
  data=data.frame(resp,period,trt,gene,sub,cov1)
  
  ## using smsn.lmm function of skewlmm package
  fm1 <- smsn.lmm(data,formFixed=(resp)~ period+ trt+gene+cov1,
                  groupVar="sub")
  summary(fm1)
  estse_sm=cbind(c( fm1$estimates$beta,fm1$estimates$sigma2,
                    fm1$estimates$D,fm1$estimates$lambda),fm1$std.error)
  festi2[sim,]=c(fm1$estimates$beta,fm1$estimates$sigma2,
                 fm1$estimates$D,fm1$estimates$lambda,fm1$std.error)
  
  ## for mahalanobis distance
  bett=fm1$estimates$beta
  sigee=fm1$estimates$sigma2
  sigss=as.numeric(fm1$estimates$D)
  si=sigee*diag(p*m)
  sigmai=si+sigss*crossprod(t(zij))
  maha_dist_smsn[sim,]=unlist(lapply(1:s, function(j) t(sapply(1:no[j], function(i) (t(yij[[j]][i,]-xij[[j]][[i]]%*%bett)%*%solve(sigmai)%*%
                                                                                       (yij[[j]][i,]-xij[[j]][[i]]%*%bett))))))
  
  ## resdiuals of smsn.lmm fit
  resd_smsn[sim,]=residuals(fm1)
  
  ## using lme function 
  fm2=lme(resp~period+trt+gene+cov1,data=data,random = ~1|sub,method = "ML")
  fm2
  ## using lmer function
  model=lmer(resp~period+trt+gene+cov1+ (1|sub),REML = F)# intercept vary randomly with the subjects
  model
  s1=summary(model)
  estse_lm=cbind(c(s1$coefficients[,1],as.numeric(VarCorr(fm2)[c(2,1),1]),0),c(
    s1$coefficients[,2],as.numeric(VarCorr(fm2)[c(2,1),2]),0)) 
  festi1[sim,]=c(s1$coefficients[,1],as.numeric(VarCorr(fm2)[c(2,1),1]),
                 s1$coefficients[,2],as.numeric(VarCorr(fm2)[c(2,1),2]))
  
  ## resiudal of lme fitting
  resd_n[sim,]=residuals(model)
  
  ##mahalanobis distance
  bett=s1$coefficients[,1]
  sigee=as.numeric(VarCorr(fm2)[c(2),1])
  sigss=as.numeric(VarCorr(fm2)[c(1),1])
  si=sigee*diag(p*m)
  sigmai=si+sigss*crossprod(t(zij))
  maha_dist_n[sim,]=unlist(lapply(1:s, function(j) t(sapply(1:no[j], function(i) (t(yij[[j]][i,]-xij[[j]][[i]]%*%bett)%*%solve(sigmai)%*%
                                                                                    (yij[[j]][i,]-xij[[j]][[i]]%*%bett))))))
  
  ## print the estimates of smsn and lme function
  cbind(test,estse_sm,estse_lm)
  
  ## aic bic of lme and smsn fit usual formula
  nofr=length(test)-1;tot=p*m*sum(no)
  aic_bic[sim,1:4]=c(-2*fm2$logLik+2*(nofr),-2*fm2$logLik+log(tot)*(nofr),
                     -2*fm1$loglik+2*(nofr),-2*fm1$loglik+log(tot)*(nofr))
  
  
  ## Now for parameter estimation EM NR method 
  beth=tbet-.1;sigmah=tsigma-.06;skewparr=lam+.05
  ## final em algorithm
  pp=parameter_estimates(beth,sigmah,skewparr)
  v0=1
  repeat{
    v0=v0+1
    bethh=pp[1:(p+m+t-2+1)]; sigmahh=pp[(p+t+m-1+1):(p+m+t+1)];skewparh=pp[(p+m+t+1+1)]
    pp=parameter_estimates(bethh,sigmahh,skewparh)
    print(cbind(pp[1:length(test)],test))
    if(max(abs(pp))>50|max(abs(pp[1:length(test)]-c(bethh,sigmahh,skewparh)))<5*10**-3) break
  }
  cat("convergence at iteration no",v0)
  
  ### for reporting in tables purposes, convert the intercept into cp
  #so we can compare it with LMM results
  la=pp[(length(tbet)+3)]
  delbb=la*as.vector((sqrt(1+la^2))^(-1))
  betcp=pp[1]+(sqrt(2/pi)*zij*sqrt(pp[(length(tbet)+2)])*delbb)[1,1]
  festi[sim,]=c(betcp,pp[-1])
  
  ## print the estimate and std of three models smsn, lm, our
  estse_our=cbind(c(betcp,pp[2:length(test)]),pp[(length(test)+1):(2*length(test))])
  print(cbind(test,estse_lm,estse_sm,estse_our))
  
  # calculate aic using log likelihood as q fn 
  bett=c(pp[1:length(tbet)]);nofr=length(test)-1
  sigee=pp[length(tbet)+1];sigss=pp[length(tbet)+2];lamm=pp[length(tbet)+3]
  si=sigee*diag(p*m)
  sigmai=si+sigss*crossprod(t(zij))
  delb=lamm/(sqrt(1+lamm^2))
  delbar=sqrt(sigss)*delb;di=zij*delbar
  bigsi=sigmai-crossprod(t(di))# same as vij
  taui=as.numeric(1/(1+t(di)%*%solve(bigsi)%*%di))
  etai=lapply(1:s, function(j) sapply(1:no[j], function(i)
    taui*t(di)%*%solve(bigsi)%*%(yij[[j]][i,]-xij[[j]][[i]]%*%bett)))
  tihat=lapply(1:s, function(j) sapply(1:no[j], function(i)
    etai[[j]][i]+sqrt(taui)*dnorm(etai[[j]][i]/sqrt(taui))/pnorm(etai[[j]][i]/sqrt(taui))))
  tisqhat=lapply(1:s, function(j) sapply(1:no[j], function(i)
    etai[[j]][i]^2+taui+etai[[j]][i]*sqrt(taui)*dnorm(etai[[j]][i]/sqrt(taui))/pnorm(etai[[j]][i]/sqrt(taui))))
  
  qfn= -.5*(sum(no)*(p*m*log(2*pi)+log(det(bigsi))-(log(2)-log(pi)))
            +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i)
              t(yij[[j]][i,]-xij[[j]][[i]]%*%bett)%*%solve(bigsi)%*%(yij[[j]][i,]-xij[[j]][[i]]%*%bett-2*tihat[[j]][i]*di)  )))) +
              Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) 
                (tisqhat[[j]][i]+t(di)%*%solve(bigsi)%*%di*tisqhat[[j]][i]))))))
  
  aic_bic[sim,5:6]=c(-2*qfn+2*nofr,-2*qfn+log(tot)*nofr)
  
   ### residual vector based on r_ij
  sn_condresid=lapply(1:s, function(j) (sapply(1:no[j], function(i) (yij[[j]][i,]-xij[[j]][[i]]%*%bett-di*tihat[[j]][i]))))
  resd_o[sim,]=unlist(sn_condresid)
  
  ## mahalanobis distance
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
festii=festi[(1:sim)[-rmm][1:ns],];festii1=festi1[(1:sim)[-rmm][1:ns],]
festii2=festi2[(1:sim)[-rmm][1:ns],]
maha_dist_nn=colMeans(maha_dist_n[(1:sim)[-rmm][1:ns],])
maha_dist_smsnn=colMeans(maha_dist_smsn[(1:sim)[-rmm][1:ns],])
maha_dist_oo=colMeans(maha_dist_o[(1:sim)[-rmm][1:ns],])
resd_nn=colMeans(resd_n[(1:sim)[-rmm][1:ns],])
resd_oo=colMeans(resd_o[(1:sim)[-rmm][1:ns],])
resd_smsnn=colMeans(resd_smsn[(1:sim)[-rmm][1:ns],])

### qqplot of mahalanobis distance
dataa=data.frame(maha_dist_nn,maha_dist_oo,maha_dist_smsnn)
colnames(dataa)=c("normal","skew error","smsn")
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
std_resd_smsn=scale(unlist(resd_smsnn))
df11=data.frame(std_resd_n,std_resd_e,std_resd_smsn)

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
fout=cbind(festii,festii1,festii2,aicc)
write.csv(fout,"D:\ ranskew out_30.csv")

## number of times our model gets selected based on AIC BIC and HQ values
length(which(aicc[,1]>(aicc[,5]-2*nofr)))/nrow(aicc)*100## AIC glosup formula


## for printing the estimate table
absbias_o= abs(rowMeans(sapply(1:nrow(festii1), function(i) abs(festii[i,1:length(test)]-test))))
absbias_n= abs(rowMeans(sapply(1:nrow(festii1), function(i) abs(festii1[i,1:length(test)]-test))))
absbias_smsn= abs(rowMeans(sapply(1:nrow(festii2), function(i) abs(festii2[i,1:length(test)]-test))))

snn=cbind(c(test),colMeans(festii[,1:length(test)]),
          colMeans(festii[,-c(1:length(test))]),
          absbias_o)
nn=cbind(colMeans(festii1[,1:(length(test)-1)]),c(test[1:(length(test)-1)]),
         colMeans(festii1[,-c(1:(length(test)-1))]),
         absbias_n[-length(test)])
smsn=cbind(colMeans(festii2[,1:length(test)]),
           colMeans(festii2[,-c(1:length(test))]),
           absbias_smsn)

ll1=cbind(snn,rbind(nn[,-2],0),smsn)
ll1
write.csv(ll1,"D:\ randomskew1_30.csv")
Sys.time()-strt

#####**** make the empirical densities of ml estimates*****
no_par=12
sim100=festii[1:ns,c(1:no_par)]

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



