

rm(list=ls())

library(survAUC)

source("RobSBoosting_fun.R") #/home2/mw846/2016tensor




model='linear'
p=1000
q=10
G_data_class='Gene'

rho=0.3
corr_str='AR'
cc_rate=20

error_type=2


if (model=='linear'){
  n=150
} else {
  n=250
}
m=100


sigma_GE=corr_setting(p,q,corr_str,rho)

NN=500

print(paste0('Set:',' Genetic type=',G_data_class,' model=',model,' n=',n))


############ training data ##############

data<-simulated_data(n,q,p,G_data_class,model,sigma_GE$sigmaG,sigma_GE$sigmaE,cc_rate,error_type=error_type)
  
alpha_true=data$para_true$alpha
beta_true=data$para_true$beta
  
  
G=data$G
E=data$E
Y=data$Y
W=data$W
E_type=data$E_type
  
  

  
  
############ testing data ##############3
  
data_test<-simulated_data(m,q,p,G_data_class,model,sigma_GE$sigmaG,sigma_GE$sigmaE,cc_rate,alpha_true,beta_true,error_type=error_type)
  
  
G_test=data_test$G
E_test=data_test$E
Y_test=data_test$Y
W_test=data_test$W

  

  
knots=list()
Boundary.knots=matrix(0,p+q,2)
for (i in 1:q){   
  knots[[i]]=c(0,1)
  Boundary.knots[i,]=c(0,1)
}


  
para_temp=c(alpha_true,beta_true)
linear_id=data$linear_id==1

  
############ KM weight for AFT model ##
  
if (model=='AFT'){
  Surv.rsp <- Surv(time=Y_cc[,1],event=Y_cc[,2])
  Surv.rsp.new <- Surv(time=Y_test[,1],event=Y_test[,2])
  w=kmw(Y[,1],Y[,2])
  delta=Y[,2]
  Y=Y[,1]
  } else {
    w=NULL
    delta=matrix(1,n,1)
    }
  



fit=boosting.robust.spline.hier(E[delta==1,],G[delta==1,],Y[delta==1],NN,num.knots=2,knots=knots,Boundary.knots=Boundary.knots,degree=2,'Robust',E_type,'nonlinear',kmweight=w[delta==1])
temp=variable_return(fit,p,q)
para_hat<-temp$para_id
unique_set=temp$unique_set
y_predict=predict_boosting(E_test,G_test,unique_set,fit)
      


temp=GetFPTP(para_temp[!linear_id],para_hat[!linear_id])
TP.nonlinear=temp$TP
      
      
temp=GetFPTP(para_temp[linear_id],para_hat[linear_id])
TP.linear=temp$TP
      
temp=GetFPTP(para_temp,para_hat)
FP=temp$FP
      
if (model=='linear'){
  PMSE=median((y_predict-Y_test)^2) 
  } else {
    PMSE<- UnoC(Surv.rsp, Surv.rsp.new, -y_predict)
  }


estimation_results=estimation_return(unique_set,fit,p,q)

estimation_results=estimation_results$estimation_results


linear_e=estimation_results[linear_id]
linear_t=para_temp[linear_id]

temp=matrix(0,sum(linear_id),1)
temp2=sapply(linear_e,'length')
temp3=which(temp2!=0)

for (ii in 1:length(temp3)){
  temp[temp3[ii]]=linear_e[[temp3[ii]]]
}

EMSE=sum((temp-linear_t)^2)/sum(linear_t!=0)



nonlinear_e=estimation_results[!linear_id]
nonlinear_t=para_temp[!linear_id]

temp=matrix(0,sum(!linear_id),1)
temp2=sapply(nonlinear_e,'length')
temp3=union(which(temp2!=0),which(nonlinear_t!=0))

xx_dot=seq(from=0,to=1,by=0.0001)

for (ii in 1:length(temp3)){
  
  bs_predict=nonlinear_e[[temp3[ii]]]
  
  if (is.null(bs_predict)){
    bs_predict=0
  }
  
  if (temp3[ii]==1){
    bs_true=6*sin(xx_dot*pi*2)-0.06
  } else if (temp3[ii]==2){
    bs_true=6*exp(2*xx_dot-1)-7.05
  } else if ((temp3[ii]==12) | (temp3[ii]==23) | (temp3[ii]==34) | (temp3[ii]==45)){
    bs_true=-2*xx_dot*(1+xx_dot)+2
  } else if ((temp3[ii]==13) | (temp3[ii]==24) | (temp3[ii]==35) | (temp3[ii]==46)){
    bs_true=-4*(xx_dot)^3+1
  } else {
    bs_true=0
  }
  
  temp[temp3[ii]]=mean((bs_true-bs_predict)^2)
}

EMISE=sum(temp)/sum(nonlinear_t!=0)


print(paste0(' TP.linear=',TP.linear,' TP.nonlinear=',TP.nonlinear,' FP=',FP,' PMSE=',round(PMSE,2),' EMSE=',round(EMSE,2),' EMISE=',round(EMISE,2)),quote = F)


vv=1

if (vv==1){
  bs_true=6*sin(xx_dot*pi*2)-0.06
} else if (vv==2){
  bs_true=6*exp(2*xx_dot-1)-7.05
} else if ((vv==12) | (vv==23) | (vv==34) | (vv==45)){
  bs_true=-2*xx_dot*(1+xx_dot)+2
} else if ((vv==13) | (vv==24) | (vv==35) | (vv==46)){
  bs_true=-4*(xx_dot)^3+1
} 


if (vv==1){
  plot(xx_dot,bs_true,type='l',lty=1,xlab='z',ylab='g(z)',ylim=c(-8,8),lwd=1.8,cex.lab=1.2,cex.axis=1.2,main='(a)')
} else if (vv==2){
  plot(xx_dot,bs_true,type='l',lty=1,xlab='z',ylab='g(z)',ylim=c(-8,8),lwd=1.8,cex.lab=1.2,cex.axis=1.2,main='(b)')
} else if ((vv==12) | (vv==23) | (vv==34) | (vv==45)){
  plot(xx_dot,bs_true,type='l',lty=1,xlab='z',ylab='g(z)',ylim=c(-3,3),lwd=1.8,cex.lab=1.2,cex.axis=1.2,main='(c)')
} else if ((vv==13) | (vv==24) | (vv==35) | (vv==46)){
  plot(xx_dot,bs_true,type='l',lty=1,xlab='z',ylab='g(z)',ylim=c(-3,3),lwd=1.8,cex.lab=1.2,cex.axis=1.2,main='(d)')
} 


lines(xx_dot,estimation_results[[vv]],col="red",lty=2,lwd=1.8)




 