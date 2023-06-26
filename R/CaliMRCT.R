#' CaliMRCT
#'
#' Calibration Weighting Estimation of Treatment Effect Heterogeneity for MRCTs
#'
#' CaliMRCT is a function to assess treatment effect heterogeneity for multi-regional randomized clinical trials.
#'
#' @param R Region indicator
#' @param X Matrix of m baseline covariates (after standardization)
#' @param A Treatment indicator
#' @param Y Observed outcome
#' @param Nboot Number of bootstrap samples used (default is 300)
#' @param seed Pseudorandom number for reproducible results (default is 2022)
#'
#' @return A list containing
#' \itemize{
#'    \item naive: a matrix containing estimated ATEs with bootstrap SEs and confidence intervals by Naive method
#'    \item ipw: a matrix containing estimated ATEs with bootstrap SEs and confidence intervals by IPW method
#'    \item cw: a matrix containing estimated ATEs with bootstrap SEs and confidence intervals by CW method
#'    \item cwe: a matrix containing estimated ATEs with bootstrap SEs and confidence intervals by CWE method
#'    \item tab: a table including results by all four methods
#' }
#' @import nleqslv
#' @export
#'
#' @author Haotian Zhuang, Xiaofei Wang <haotian.zhuang@@duke.edu>
#' @examples CaliMRCT(R = simMRCT$R, X = cbind(simMRCT$X1,simMRCT$X2), A = simMRCT$A, Y = simMRCT$Y, Nboot = 300, seed = 2022)

CaliMRCT<-function(R,X,A,Y,Nboot=300,seed=2022){

  set.seed(seed = seed)

  X1<-X[R==1,]
  X2<-X[R==2,]
  A1<-A[R==1]
  A2<-A[R==2]
  Y1<-Y[R==1]
  Y2<-Y[R==2]

  beta1_unadj=mean(Y1[A1==1])-mean(Y1[A1==0])
  beta2_unadj=mean(Y2[A2==1])-mean(Y2[A2==0])

  w1<-IPW(X1 = X1,X2 = X2)$W1
  w2<-IPW(X1 = X1,X2 = X2)$W2
  beta1_ipw=sum(w1*(A1*Y1))/sum(w1[A1==1])-sum(w1*((1-A1)*Y1))/sum(w1[A1==0])
  beta2_ipw=sum(w2*(A2*Y2))/sum(w2[A2==1])-sum(w2*((1-A2)*Y2))/sum(w2[A2==0])

  p=weight(X = X1,gtilde = apply(cbind(rbind(X1,X2),rbind(X1,X2)^2),2,mean))
  q=weight(X = X2,gtilde = apply(cbind(rbind(X1,X2),rbind(X1,X2)^2),2,mean))

  beta1_cw=sum(p*(A1*Y1-(1-A1)*Y1))/0.5
  beta2_cw=sum(q*(A2*Y2-(1-A2)*Y2))/0.5
  beta1_cwe=sum(p*(A1*Y1))/sum(p[A1==1])-sum(p*((1-A1)*Y1))/sum(p[A1==0])
  beta2_cwe=sum(q*(A2*Y2))/sum(q[A2==1])-sum(q*((1-A2)*Y2))/sum(q[A2==0])

  unadjDiff<-beta2_unadj-beta1_unadj
  ipwDiff<-beta2_ipw-beta1_ipw
  cwDiff<-beta2_cw-beta1_cw
  cweDiff<-beta2_cwe-beta1_cwe

  unadjPool<-length(Y1)/(length(Y1)+length(Y2))*beta1_unadj+length(Y2)/(length(Y1)+length(Y2))*beta2_unadj
  ipwPool<-length(Y1)/(length(Y1)+length(Y2))*beta1_ipw+length(Y2)/(length(Y1)+length(Y2))*beta2_ipw
  cwPool<-length(Y1)/(length(Y1)+length(Y2))*beta1_cw+length(Y2)/(length(Y1)+length(Y2))*beta2_cw
  cwePool<-length(Y1)/(length(Y1)+length(Y2))*beta1_cwe+length(Y2)/(length(Y1)+length(Y2))*beta2_cwe

  datb<-list(X1=X1,X2=X2,A1=A1,A2=A2,Y1=Y1,Y2=Y2)
  tebt<-teboot(Nboot = Nboot,dat = datb)

  all<-c(beta1_unadj,beta2_unadj,unadjDiff,unadjPool,
         beta1_ipw,beta2_ipw,ipwDiff,ipwPool,
         beta1_cw,beta2_cw,cwDiff,cwPool,
         beta1_cwe,beta2_cwe,cweDiff,cwePool)

  method=c(rep("Naive",4),rep("IPW",4),rep("CW",4),rep("CWE",4))
  ATE=c(rep(c("Region 1","Region 2","Difference","Pool"),4))
  tab=cbind(method,ATE,round(cbind(all,tebt$bsd,tebt$cil,tebt$ciu),3))
  colnames(tab)<-c("Method","ATE","Est","SE","95% LB","95% UB")
  rownames(tab)<-NULL

  results<-cbind(all,tebt$bsd,tebt$cil,tebt$ciu)
  rownames(results)<-ATE
  colnames(results)<-c("Est","SE","95% LB","95% UB")
  naive<-results[1:4,]
  ipw<-results[5:8,]
  cw<-results[9:12,]
  cwe<-results[13:16,]

  return(list(naive=naive,ipw=ipw,cw=cw,cwe=cwe,tab=tab))

}


weight<-function(X,gtilde){
  n<-nrow(X)
  k<-length(gtilde)

  fun<-function(lambda){
    f<-sapply(1:k, function(i) {
      sum(exp(cbind(X,X^2)%*%lambda[1:k])*(cbind(X,X^2)[,i]-gtilde[i]))
    })
  }

  sol<-nleqslv(rep(0,k),fun)
  lambda<-sol$x

  p<-exp(cbind(X,X^2)%*%lambda[1:k])/sum(exp(cbind(X,X^2)%*%lambda[1:k]))
}

# IPW method
expit <- function(x) {return (exp(x)/(1+exp(x)))}

IPW<-function(X1,X2){
  n1<-nrow(X1)
  n2<-nrow(X2)
  r1<-c(rep(1,n1),rep(0,n2))
  r2<-c(rep(1,n2),rep(0,n1))

  s1<-glm(r1 ~ rbind(X1,X2),family = binomial("logit"))
  w1<-1/as.numeric(expit(cbind(1,X1)%*%s1$coefficients))

  s2<-glm(r2 ~ rbind(X2,X1),family = binomial("logit"))
  w2<-1/as.numeric(expit(cbind(1,X2)%*%s2$coefficients))

  return(list(W1=w1,W2=w2))
}

# Bootstrap estimation
teboot<-function(Nboot=300,dat){
  B<-matrix(NA,Nboot,16)
  N1<-nrow(dat$X1)
  N2<-nrow(dat$X2)
  for (kkk in 1:Nboot) {
    B1=sample(1:N1,N1,replace = T)
    B2=sample(1:N2,N2,replace = T)
    X1<-as.matrix(dat$X1[B1,])
    X2<-as.matrix(dat$X2[B2,])
    A1<-dat$A1[B1]
    A2<-dat$A2[B2]
    Y1<-dat$Y1[B1]
    Y2<-dat$Y2[B2]

    p=weight(X = X1,gtilde = apply(cbind(rbind(X1,X2),rbind(X1,X2)^2),2,mean))
    q=weight(X = X2,gtilde = apply(cbind(rbind(X1,X2),rbind(X1,X2)^2),2,mean))

    B[kkk,1]<-(mean(Y1[A1==1])-mean(Y1[A1==0]))
    B[kkk,2]<-(mean(Y2[A2==1])-mean(Y2[A2==0]))
    B[kkk,3]<-B[kkk,2]-B[kkk,1]
    B[kkk,4]<-N1/(N1+N2)*B[kkk,1]+N2/(N1+N2)*B[kkk,2]

    w1<-IPW(X1 = X1,X2 = X2)$W1
    w2<-IPW(X1 = X1,X2 = X2)$W2
    B[kkk,5]=sum(w1*(A1*Y1))/sum(w1[A1==1])-sum(w1*((1-A1)*Y1))/sum(w1[A1==0])
    B[kkk,6]=sum(w2*(A2*Y2))/sum(w2[A2==1])-sum(w2*((1-A2)*Y2))/sum(w2[A2==0])
    B[kkk,7]=B[kkk,6]-B[kkk,5]
    B[kkk,8]<-N1/(N1+N2)*B[kkk,5]+N2/(N1+N2)*B[kkk,6]

    B[kkk,9]<-sum(p*(A1*Y1-(1-A1)*Y1))/0.5
    B[kkk,10]<-sum(q*(A2*Y2-(1-A2)*Y2))/0.5
    B[kkk,11]<-B[kkk,10]-B[kkk,9]
    B[kkk,12]<-N1/(N1+N2)*B[kkk,9]+N2/(N1+N2)*B[kkk,10]

    B[kkk,13]<-(sum(p*(A1*Y1))/sum(p[A1==1])-sum(p*((1-A1)*Y1))/sum(p[A1==0]))
    B[kkk,14]<-(sum(q*(A2*Y2))/sum(q[A2==1])-sum(q*((1-A2)*Y2))/sum(q[A2==0]))
    B[kkk,15]<-B[kkk,14]-B[kkk,13]
    B[kkk,16]<-N1/(N1+N2)*B[kkk,13]+N2/(N1+N2)*B[kkk,14]
  }
  bsd<-apply(B,2,sd)
  return(list(bsd=bsd,cil=apply(B,2,quantile,probs=0.025),ciu=apply(B,2,quantile,probs=0.975)))
}
