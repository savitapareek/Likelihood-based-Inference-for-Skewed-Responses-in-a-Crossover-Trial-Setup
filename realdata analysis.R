## cleaning the enviornment
rm(list=ls())
gc()

## running time
strt=Sys.time()

## install requisite libraries
library(psych)# for trace of matrix
library(MASS)# for mvtnorm
library(lme4) # for lme4 function
library(magic)## for adiag
library(moments)## for skewness
library(nlme) # for lme function 
library(skewlmm)## for smsn.lmm function
library(matlib)# ginv R
library(EnvStats) # for qqPlot
library(ggplot2)
library(ggpubr)

#######**********error skewnormal parameter estimate function%%%%******&&&&&&####$$$$
parameter_estimates=function(beth,sigmah,skewparr){
  sigmaeh=sigmah[1]; sigmash=sigmah[2];skewparh=c(skewparr,rep(0,((p*m)-1)))
  delij=skewparh*as.vector((sqrt(1+crossprod(skewparh)))^(-1))
  muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]]%*%beth))
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
  bethh=(solve(Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(xij[[j]])%*%solve(vij)%*%xij[[j]])))))
         %*%Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(xij[[j]])%*%solve(vij)%*%(yij[[j]][i,]-dij*T01[[j]][i]))))))
  
  ## for information matrix
  sbet= Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]])%*%solve(vij)%*%xij[[j]]))))
  sebet=sqrt(diag(solve(-sbet)))
  vo=1
  repeat{
    vo=vo+1
    sigp=c(sigmaeh,sigmash,skewparr);skewparh=c(skewparr,rep(0,(p*m-1)))
    delij=skewparh*as.vector((sqrt(1+crossprod(skewparh)))^(-1))
    muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]]%*%bethh))
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
    
    ## derivative of dij with respect to sigma_e^2, sigma_s^2,  \lambda
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
    #if(R(H)<nopar) sigmahh=c(sigp)-Ginv(H)%*%dQ_esvrl
    # if(R(H)==nopar) 
    sigmahh=c(sigp)-solve(H)%*%dQ_esvrl
    #print(cbind(sigmahh,c(isigma,skewparr)))
    sigmaeh=sigmahh[1];sigmash=sigmahh[2];skewparr=sigmahh[3]
    if(max(abs(sigmahh-sigp))<5*10^-3) break #max(abs(sigmahh))>50|vo>500|
  }
  #print(cbind(sigmahh,c(tsigma,lam)))
  #sevar= sqrt(diag(solve(-H)))
  
  ## for standard error as diagonal elements of inverse ofInformation matrix
  dQbet_esvrl=sapply(1:nopar, function(k)  -.5*(t(Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i)  t(yij[[j]][i,]-muij[[j]][[i]])%*%(dvijinv[[k]])%*%(-xij[[j]]))))))
                                                +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]])%*%dvijinv[[k]]%*%(yij[[j]][i,]-muij[[j]][[i]]-2*dij*T01[[j]][i])) )))
                                                +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]])%*%solve(vij)%*%(-2*ddij[[k]]*T01[[j]][i])) ))))
  )
  
  I= matrix(nrow=length(beth)+3,ncol=length(beth)+3)
  pm=(p+t+m-2)
  I[1:pm,1:pm]=sbet
  for (i in 1:nopar) {
    I[1:pm,(pm+i)]=I[(pm+i),1:pm]=dQbet_esvrl[,i]
  }
  I[((pm+1):(pm+nopar)),((pm+1):(pm+nopar))]=H
  ## for standard error as diagonal elements of inverse of Information matrix
  # if(R(I)<length(beth)+3)  see=sqrt(abs(diag(Ginv(-I))))
  #if(R(I)==length(beth)+3) 
  see=sqrt(abs(diag(solve(-I))))
  return(c(bethh,sigmahh,see)) 
  
}

###********#######%%#$$$$ for random effect skewnorm parameter estimate function ****$$$$%%%%#####
parameter_estimates1=function(beth,sigmah,skewparr){
  sigmaeh=sigmah[1]; sigmash=sigmah[2]
  delb=skewparr*as.vector((sqrt(1+skewparr^2))^(-1))
  muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]]%*%beth))
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
  bethh=(solve(Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(xij[[j]])%*%solve(vij)%*%xij[[j]])))))
         %*%Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) t(xij[[j]])%*%solve(vij)%*%(yij[[j]][i,]-dij*T01[[j]][i]))))))
  
  ## for information matrix of beta
  sbet= Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]])%*%solve(vij)%*%xij[[j]]))))
  
  ## use NR method for estimating variance and skewness parameters
  repeat{
    sigp=c(sigmaeh,sigmash,skewparr)
    delb=skewparr*as.vector((sqrt(1+skewparr^2))^(-1))
    muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]]%*%bethh))
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
    
    ## derivative of Q function with respect to sigma_e^2, sigma_s^2, \lambda
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
    
    ## second derivative of solve(vij) with respect to sigma_e^2, sigma_s^2,  \lambda
    ddvijinv=  lapply(1:nopar, function(u) lapply(u:nopar, function(u1)
      (-2*solve(vij)%*%dvij[[u]]%*%dvijinv[[u1]]
       -solve(vij)%*%ddvij[[u]][[u1-u+1]]%*%solve(vij))))
    
    ## second derivative of Q function with respect to sigma_e^2, sigma_s^2,  \lambda
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
    #print(cbind(sigmahh,c(tsigma,lam)))
    sigmaeh=sigmahh[1];sigmash=sigmahh[2];skewparr=sigmahh[3]
    if(max(abs(sigmahh-sigp))<5*10^-3) break #max(abs(sigmahh))>50|
  }
  
  dQbet_esvrl=sapply(1:nopar, function(k)  -.5*(t(Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i)  t(yij[[j]][i,]-muij[[j]][[i]])%*%(dvijinv[[k]])%*%(-xij[[j]]))))))
                                                +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]])%*%dvijinv[[k]]%*%(yij[[j]][i,]-muij[[j]][[i]]-2*dij*T01[[j]][i])) )))
                                                +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) -t(xij[[j]])%*%solve(vij)%*%(-2*ddij[[k]]*T01[[j]][i])) ))))
  )
  
  I= matrix(nrow=length(beth)+3,ncol=length(beth)+3)
  pm=(p+t+m-2)
  I[1:pm,1:pm]=sbet
  for (i in 1:nopar) {
    I[1:pm,(pm+i)]=I[(pm+i),1:pm]=dQbet_esvrl[,i]
  }
  I[((pm+1):(pm+nopar)),((pm+1):(pm+nopar))]=H
  ## for standard error as diagonal elements of inverse of Information matrix
  # if(R(I)<length(beth)+3)  see=sqrt(abs(diag(Ginv(-I))))
  # if(R(I)==length(beth)+3)  
  see=sqrt(abs(diag(solve(-I))))
  return(c(bethh,sigmahh,see))
  
}

## now define the models specifications 
s=3#no of sequences
p=t=3#no of periods and no of treatments
m=10#no of  genes
no=array()
pm=p*m
no=c(4,4,4)
treatment1=list() ## we assume placebo as treatment 1, 15mg as treatment 2 and 25mg as treatment 3
treatment1[[1]]=c(2,1,3)
treatment1[[2]]=c(3,2,1)
treatment1[[3]]=c(1,3,2)

######################### treatment matrix for xij ###
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


############## xij matrix
x1k=list()
for (i in 1:s) {
  x1k[[i]]=cbind(rep(1,p*m),rbind(matrix(data=0,nrow = m,ncol = p-1),kronecker(diag(p-1),rep(1,m))),Trt[[i]],kronecker(rep(1,p),rbind(rep(0,m-1),diag(m-1))))#,carry[[i]]
}

#making zij matrix
z1k=cbind(kronecker(rep(1,p),rep(1,m)))

###Read the gene expression real data, #anova 3 way all seq excel file
data11=read.csv(file.choose(),header = T)

### construct the subject and sequence categorical variable
sub=rep(paste("sub",1:sum(no),sep=""),each=p*m)
seq=rep(c("seq1","seq2","seq3"),each=p*m*no[1])

### data frame to be used for lme and smsn.lmm fitting
data=data.frame(data11,sub,seq)
attach(data)
head(data)

###*** fit the gene data using skewlmm package with smsn.lmm function****************
fm1 <- smsn.lmm(data,formFixed=(Response_allseq)~ Period+ Treatment+Gene,
                groupVar="sub")
summary(fm1)
plot(fm1)
healy.plot(fm1)
s1=summary(fm1)
coef_smsnnorm=c(s1$tableFixed[c(1:5,7:14,6)],(s1$varFixed),(s1$varRandom),
                fm1$estimates$lambda,s1$criteria[2:3],sqrt(mean(residuals(fm1)^2)))

###**********fit the gene data using lmer function************
model=lmer((Response_allseq)~ Period+ Treatment+Gene+ (1|sub),REML = F)# intercept vary randomly with the subjects
model
s1=summary(model)
plot(model)
##  compute the root mean square error RMSE
rmse_lme=sqrt(mean(residuals(model)^2))#or sqrt(var(residuals(model)))
## fixed and random effects
re_dat = as.data.frame(VarCorr(model))
int_vcov = re_dat[1,'vcov']
resid_vcov = re_dat[2,'vcov']
bet_org=c(s1$coefficients[,1][-6],s1$coefficients[,1][6])
betse=c(s1$coefficients[,2][-6],s1$coefficients[,2][6])
options(scipen=10)
tot=p*m*sum(no);nofr=p+t+m-2+2
nse=c(betse,re_dat[1,'sdcor'],re_dat[2,'sdcor'],rep(0,4))
coef_lme=cbind(c(bet_org,s1$sigma**2,int_vcov,0,extractAIC(model)[2],
                 -2*s1$logLik+nofr*log(tot),rmse_lme),nse)
cbind(coef_lme,coef_smsnnorm)
#write.csv(cbind(coef_lme,coef_smsnnorm),"D:\ out.csv")

######******value of y and initial parameters for fitting the proposed algorithm
y=list()
y[[1]]=data11[(1:120),1]
y[[2]]=data11[(121:240),1]
y[[3]]=data11[(241:360),1]

isigmae=resid_vcov
isigmas=int_vcov
ibet=bet_org
skewparr=.001
isigma=c(isigmae,isigmas)

## xij yij and zij are
xij=x1k;zij=z1k
yij=lapply(1:s, function(i) matrix(y[[i]],ncol=30,byrow=T))

### ## compute mahalanobis distance ddi for lme fit
vij=isigmae*diag(p*m)+isigmas*crossprod(t(zij))
ddin=unlist(lapply(1:s, function(j) sapply(1:no[j], function(i) (t(yij[[j]][i,]-xij[[j]]%*%bet_org)%*%solve(vij)%*%
                                                                   (yij[[j]][i,]-xij[[j]]%*%bet_org)))))

##************** fit the random error to be skewed#####******%%%!!!%%%%%@@@@******************
pp=parameter_estimates(ibet,isigma,skewparr)
pp[1:17]
v0=1
repeat{
  v0=v0+1
  beth=pp[1:(p+m+t-2)]; sigmah=pp[(p+t+m-1):(p+m+t)];skewparh=pp[(p+m+t+1)]
  pp=parameter_estimates(beth,sigmah,skewparh)
  print(pp)
  if(max(abs(pp[1:17]-c(beth,sigmah,skewparh)))<5*10**-3) break
}
v0
coef_errskew=cbind(pp[c(1:17)],pp[-c(1:17)])

### calculate Tij and otehr quantites of Q-fn
skewparre=pp[17];bethhe=pp[1:14];sigmashe=pp[16];sigmaehe=pp[15]
skewparh=c(skewparre,rep(0,((p*m)-1)))
delij=skewparh*as.vector((sqrt(1+crossprod(skewparh)))^(-1))
muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]]%*%bethhe))
dij=sqrt(sigmaehe)*(delij)
R=(diag(p*m)-crossprod(t(delij)))
vij=sigmaehe*R+sigmashe*crossprod(t(zij))
tausq=as.numeric(1/(1+(t(dij)%*%solve(vij)%*%dij)))
eta=lapply(1:s, function(j) sapply(1:no[j], function(i) t(dij)%*%solve(vij)%*%(yij[[j]][i,]-muij[[j]][,i])*tausq))
T01=lapply(1:s, function(j) sapply(1:no[j], function(i) eta[[j]][i]+dnorm(eta[[j]][i]/sqrt(tausq))/pnorm(eta[[j]][i]/sqrt(tausq))
                                   *sqrt(tausq)))

## conditional residuals using model eq having t
resd_e=lapply(1:s, function(j) (sapply(1:no[j], function(i) (yij[[j]][i,]-muij[[j]][,i]-dij*T01[[j]][i]))))
rm_er=sqrt(mean((unlist(resd_e))^2))
shapiro.test(unlist(resd_e))
x2=unlist(resd_e)[-262]
ks.test(x2, "pnorm", mean=mean(x2), sd=sd(x2))
### QQ plot of std residual with ggplot removing the 262nd observation
std_resd_e=scale(unlist(resd_e))[-262]
df1=data.frame(std_resd_e)
fig1=ggplot(df1, aes(sample = std_resd_e)) + stat_qq(shape=1) +stat_qq_line()+
  theme_classic()+
  xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  theme(axis.title=element_text(size=10),plot.title = element_text(size=12))+
  ggtitle("Normal Q-Q plot of \n standardized residuals")
fig1

## mahalanobis distance ## compute mahalanobis distance ddi
v2ij=sigmaehe*diag(p*m)+sigmashe*crossprod(t(zij))
ddies=unlist(lapply(1:s, function(j) sapply(1:no[j], function(i) (t(yij[[j]][i,]-xij[[j]]%*%bethhe)%*%solve(v2ij)%*%
                                                                    (yij[[j]][i,]-xij[[j]]%*%bethhe)))))

## Q-Q plot for Chi^2 data against true theoretical data after removing 9th subject
dataa=data.frame(ddies[-9])
fig2 <- ggplot(dataa, aes(sample = ddies..9.)) + 
  stat_qq(distribution = stats::qchisq,dparams = p*m,shape=1) + 
  stat_qq_line(distribution = stats::qchisq,dparams = p*m)+ theme_classic()+
  ggtitle("Chi-square Q-Q plot of \n Mahalanobis distances")+
  theme(axis.title=element_text(size=10),plot.title = element_text(size=12))+
  xlab("Theoretical Quantiles") + ylab("Sample Quantiles")
fig2
f1=ggarrange(fig2,fig1)
f1
ks.test(dataa,pchisq,df=p*m)
##for aic and bic formula
bett=pp[1:length(ibet)];nofr=length(ibet)+2
sigee=pp[length(ibet)+1];sigss=pp[length(ibet)+2];lamm=pp[length(ibet)+3]
lambr=c(lamm, rep(0,(pm-1)))
si=sigee*diag(p*m)#;bett=pp[1:length(tbet)]
sigmai=si+sigss*crossprod(t(zij))
delb=lambr*as.vector(sqrt((1+crossprod(lambr)))^(-1))
delbar=sqrt(sigee)*delb;di=delbar
bigsi=sigmai-crossprod(t(di))# same as vij
taui=as.numeric(1/(1+t(di)%*%solve(bigsi)%*%di))
etai=lapply(1:s, function(j) sapply(1:no[j], function(i)
  taui*t(di)%*%solve(bigsi)%*%(yij[[j]][i,]-xij[[j]]%*%bett)))
tihat=lapply(1:s, function(j) sapply(1:no[j], function(i)
  etai[[j]][i]+sqrt(taui)*dnorm(etai[[j]][i]/sqrt(taui))/pnorm(etai[[j]][i]/sqrt(taui))))
tisqhat=lapply(1:s, function(j) sapply(1:no[j], function(i)
  etai[[j]][i]^2+taui+etai[[j]][i]*sqrt(taui)*dnorm(etai[[j]][i]/sqrt(taui))/pnorm(etai[[j]][i]/sqrt(taui))))

qfn= -.5*(sum(no)*(p*m*log(2*pi)+log(det(bigsi))-(log(2)-log(pi)))+Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i)
  t(yij[[j]][i,]-xij[[j]]%*%bett)%*%solve(bigsi)%*%(yij[[j]][i,]-xij[[j]]%*%bett-2*tihat[[j]][i]*di)  )))) +
    Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i)
      (tisqhat[[j]][i]+t(di)%*%solve(bigsi)%*%di*tisqhat[[j]][i]))))))

aic_er=-2*qfn+2*nofr#-qfn/(p*m*sum(no))+(nofr/(p*m*sum(no)))
bic_er=-2*qfn+log(tot)*nofr

## final output
coef_errskew=cbind(c(pp[c(1:17)],aic_er,bic_er,rm_er),c(pp[-c(1:17)],rep(0,3)))
fout=cbind(coef_lme,coef_smsnnorm,coef_errskew)

##****&&&&%%%###@@@ when b or random effect is skew normal***********&&&&&&&&&&&&&&$$$$$$$$$$$$$$$$$$$$
skewparr=.001;isigma=c(isigmae,isigmas);ibet=bet_org
pp=parameter_estimates1(ibet,isigma,skewparr)
pp[1:17]
v0=1
repeat{
  v0=v0+1
  beth=pp[1:(p+m+t-2)]; sigmah=pp[(p+t+m-1):(p+m+t)];skewparh=pp[(p+m+t+1)]
  pp=parameter_estimates1(beth,sigmah,skewparh)
  print(pp)
  if(max(abs(pp[1:17]-c(beth,sigmah,skewparh)))<5*10**-3) break
}
v0 ## algorithm converges in vo steps
coef_ran=cbind(pp[1:17],pp[-c(1:17)])

## calculate Tij and other componnets of Q-function
skewparrb=pp[17];bethhb=pp[1:14];sigmashb=pp[16];sigmaehb=pp[15]
delb=skewparrb*as.vector((sqrt(1+skewparrb^2))^(-1))
muij=lapply(1:s, function(i) sapply(1:no[i], function(j) xij[[i]]%*%bethhb))
dij=sqrt(sigmashb)*(delb)*zij
R=(1-delb^2)
vij=sigmaehb*diag(p*m)+R*sigmashb*crossprod(t(zij))
tausq=as.numeric(1/(1+(t(dij)%*%solve(vij)%*%dij)))
eta=lapply(1:s, function(j) sapply(1:no[j], function(i) t(dij)%*%solve(vij)%*%(yij[[j]][i,]-muij[[j]][,i])*tausq))
T01=lapply(1:s, function(j) sapply(1:no[j], function(i) eta[[j]][i]+dnorm(eta[[j]][i]/sqrt(tausq))/pnorm(eta[[j]][i]/sqrt(tausq))
                                   *sqrt(tausq)))
resd_b=lapply(1:s, function(j) (sapply(1:no[j], function(i) (yij[[j]][i,]-muij[[j]][,i]-dij*T01[[j]][i]))))
rm_b=sqrt(mean((unlist(resd_b))^2))

## qqplot of standardized residual after removing 262nd observation
std_resd_b=scale(unlist(resd_b))[-262]
df1=data.frame(std_resd_b)
fig3=ggplot(df1, aes(sample = std_resd_b)) + stat_qq() +stat_qq_line()+
  theme_classic()+
  xlab("Theoretical Quantiles")+ylab("Sample Quantiles")+
  theme(axis.title=element_text(size=10),plot.title = element_text(size=12))+
  ggtitle("Q-Q plot of the vector r_ij")
fig3

## for mahalanobis distance d_ij
v2ij=sigmaehb*diag(p*m)+sigmashb*crossprod(t(zij))
ddibs=unlist(lapply(1:s, function(j) sapply(1:no[j], function(i) (t(yij[[j]][i,]-xij[[j]]%*%bethhb)%*%solve(v2ij)%*%
                                                                    (yij[[j]][i,]-xij[[j]]%*%bethhb)))))

## chi square qqplot of mahalanobis distance after removing 9th subject
dataa=data.frame(ddibs[-9])
fig4 <- ggplot(dataa, aes(sample = ddibs..9.)) +
  stat_qq(distribution = stats::qchisq,dparams = p*m) + 
  stat_qq_line(distribution = stats::qchisq,dparams = p*m)+ theme_classic()+
  ggtitle("Q-Q plot of Mahalanobis\n distances")+
  theme(axis.title=element_text(size=10),plot.title = element_text(size=12))+
  xlab("Theoretical Quantiles") + ylab("Sample Quantiles")
fig4
f2=ggarrange(fig3,fig4)
f2
ggarrange(f1,f2,nrow = 2)

# calculate aic using log likelihood as q fn 
bett=c(pp[1:length(ibet)]);nofr=length(ibet)+2
sigee=pp[length(ibet)+1];sigss=pp[length(ibet)+2];lamm=pp[length(ibet)+3]
si=sigee*diag(p*m)#;bett=pp[1:length(tbet)]
sigmai=si+sigss*crossprod(t(zij))
delb=lamm/(sqrt(1+lamm^2))
delbar=sqrt(sigss)*delb;di=zij*delbar
bigsi=sigmai-crossprod(t(di))# same as vij
taui=as.numeric(1/(1+t(di)%*%solve(bigsi)%*%di))
etai=lapply(1:s, function(j) sapply(1:no[j], function(i)
  taui*t(di)%*%solve(bigsi)%*%(yij[[j]][i,]-xij[[j]]%*%bett)))
tihat=lapply(1:s, function(j) sapply(1:no[j], function(i)
  etai[[j]][i]+sqrt(taui)*dnorm(etai[[j]][i]/sqrt(taui))/pnorm(etai[[j]][i]/sqrt(taui))))
tisqhat=lapply(1:s, function(j) sapply(1:no[j], function(i)
  etai[[j]][i]^2+taui+etai[[j]][i]*sqrt(taui)*dnorm(etai[[j]][i]/sqrt(taui))/pnorm(etai[[j]][i]/sqrt(taui))))

qfn= -.5*(sum(no)*(p*m*log(2*pi)+log(det(bigsi))-(log(2)-log(pi)))
          +Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i)
            t(yij[[j]][i,]-xij[[j]]%*%bett)%*%solve(bigsi)%*%(yij[[j]][i,]-xij[[j]]%*%bett-2*tihat[[j]][i]*di)  )))) +
            Reduce("+",lapply(1:s, function(j) Reduce("+",lapply(1:no[j], function(i) 
              (tisqhat[[j]][i]+t(di)%*%solve(bigsi)%*%di*tisqhat[[j]][i]))))))

aic_ra=-2*qfn+2*nofr#(-qfn/(p*m*sum(no)))+(nofr/(p*m*sum(no)))
bic_ra=-2*qfn+log(tot)*nofr
aic_ra;bic_ra
coef_randskew=cbind(c(pp[1:17],aic_ra,bic_ra,rm_b),c(pp[-c(1:17)],rep(0,3)))

## final out put of all 4 models and store them in a csv file
fout1=cbind(coef_lme,coef_smsnnorm,coef_randskew,coef_errskew)
fout1
write.csv(fout1,"D:\ real skew.csv")

###***** plot these mahalanobis distances in one chart *******
dataa=data.frame(ddin,ddibs,ddies)[-9,]
colnames(dataa)=c("normal","skew random","skew error")
dataa

pl <- ggplot(dataa, aes(sample = dataa[,1]))+
  stat_qq(distribution = stats::qchisq,dparams = 30) + 
  stat_qq_line(distribution = stats::qchisq,dparams = 30)+
  stat_qq(aes(sample = dataa[,3]),col="blue",distribution = stats::qchisq,dparams = 30)+
  stat_qq(aes(sample = dataa[,2]),col="black",distribution = stats::qchisq,dparams = 30)+
  theme_classic()
pl

probDist <- pchisq(dataa[,1],df=30)
probDist1 <- pchisq(dataa[,2],df=30)
probDist2 <- pchisq(dataa[,3],df=30)
#create PP plot
plot(ppoints(length(dataa[,1])), sort(probDist),col="black",type="l", main = 'PP Plot', xlab = 'Observed Probability', ylab = 'Expected Probability',ylim=c(0,1))
#add diagonal line
abline(0,1)
lines(ppoints(length(dataa[,1])), sort(probDist1),col="black")#, main = 'PP Plot', xlab = 'Observed Probability', ylab = 'Expected Probability')
lines(ppoints(length(dataa[,2])), sort(probDist2),col="blue")
