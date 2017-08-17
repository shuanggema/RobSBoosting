rm(list=ls())


library(MASS)
library('survAUC')

source("RobSBoosting_fun.R") #/home2/mw846/2016tensor


model='linear'
p=1000
q=10
G_data_class='Gene'

rho=0.3
corr_str='AR'
cc_rate=20

error_type=1


if (model=='linear'){
  n=150
} else {
  n=250
}
m=100


sigma_GE=corr_setting(p,q,corr_str,rho)

error_type=1

print(paste0('Set:',' Genetic type=',G_data_class,' model=',model,' n=',n))

NN=500
tau=0.3

miss_number=1
miss_rate=20





result=list()




print(paste0('Set:',' Genetic type=',G_data_class,' model=',model,' miss_number=',miss_number,' miss_rate=',miss_rate))


  
#######completed training data############
data<-simulated_data(n,q,p,G_data_class,model,sigma_GE$sigmaG,sigma_GE$sigmaE,cc_rate,error_type=error_type)
  
alpha_true=data$para_true$alpha
beta_true=data$para_true$beta
  
  
G=data$G
E=data$E
Y=data$Y
W=data$W
E_type=data$E_type
  
  
  
  
E_cc=E
E_full=E
Y_cc=Y
  
  
############ testing data ##############3
  
data_test<-simulated_data(m,q,p,G_data_class,model,sigma_GE$sigmaG,sigma_GE$sigmaE,cc_rate,alpha_true,beta_true,error_type=error_type)
  

G_test=data_test$G
E_test=data_test$E
Y_test=data_test$Y
W_test=data_test$W

  
  
####### missing data generation ###########
  
  
  
  
  
  
miss_id=1:miss_number
nonmiss_id=setdiff(1:q,miss_id)
if (miss_number==1){
  nonmiss_id1=c(2,3,4,7,8)
  } else if (miss_number==2){
    nonmiss_id1=c(3,4,5,6,7,8,9,10)
    }
  
  
EE=cbind(E_cc[,nonmiss_id1])
p1=dim(EE)[2]
  

  
miss_each=matrix(c(1,2),1,2)
  
if (model=='linear'){
    rate=0
    con=T
    while (con){
      miss_each=matrix(0,n,miss_number)
      u=matrix(0,n,miss_number)
      for (ii in 1:miss_number){
          alpha=0.6*runif(p1,0,1)
          if (miss_number==1) {
            if (miss_rate==20){
              u[,ii]=1/(1+exp(-(EE%*%alpha+1)))  # 15: 2.5 20: 1 30:0.5
            } else if (miss_rate==40){
              u[,ii]=1/(1+exp(-(EE%*%alpha+0.5)))  # 15: 2.5 20: 1 30:0.5
            }
          } else {
            if (miss_rate==20){
              u[,ii]=1/(1+exp(-(EE%*%alpha+1)))  # 15: 2.5 20: 1 30:0.5
            } else if (miss_rate==40){
              u[,ii]=1/(1+exp(-(EE%*%alpha+0.5)))  # 15: 2.5 20: 1 30:0.5
            }
          }
        }
        
        a=matrix(runif(n*miss_number,0,1),n,miss_number)
        
        miss_each[u>a]=1# observed as 1; missing as 0
        
        miss_temp=rowSums(miss_each)
        miss=matrix(0,n,1)
        miss[miss_temp==miss_number]=1
        rate=sum(miss==0)/n# missing rate
        
        if (miss_number==1){
          con=abs(rate-miss_rate*0.01)>0.02
        } else {
          con=((abs(rate-miss_rate*0.01)>0.02) | (sum(miss_each[,1])>sum(miss_each[,2])))
        }
        
      }
      for (ii in 1:miss_number){
        E_cc[miss_each[,ii]==0,miss_id[ii]]=NA
      }

  } else {
    rate=0
    
    
    dead_id=which(Y[,2]==1)
    EE=EE[dead_id,]
    n_dead=sum(Y[,2]==1)
    
    con=T
    while (con){
      miss_each=matrix(0,n_dead,miss_number)
      u=matrix(0,n_dead,miss_number)
      for (ii in 1:miss_number){
        alpha=0.6*runif(p1,0,1)
        if (miss_number==1) {
          if (miss_rate==20){
            u[,ii]=1/(1+exp(-(EE%*%alpha+1)))  # 15: 2.5 20: 1 30:0.5
          } else if (miss_rate==40){
            u[,ii]=1/(1+exp(-(EE%*%alpha+0.5)))  # 15: 2.5 20: 1 30:0.5
          }
        } else {
          if (miss_rate==20){
            u[,ii]=1/(1+exp(-(EE%*%alpha+1)))  # 15: 2.5 20: 1 30:0.5
          } else if (miss_rate==40){
            u[,ii]=1/(1+exp(-(EE%*%alpha+0.5)))  # 15: 2.5 20: 1 30:0.5
          }
        }
      }
      
      a=matrix(runif(n_dead*miss_number,0,1),n_dead,miss_number)
      
      miss_each[u>a]=1# observed as 1; missing as 0
      
      miss_temp=rowSums(miss_each)
      miss=matrix(0,n_dead,1)
      miss[miss_temp==miss_number]=1
      rate=sum(miss==0)/n# missing rate
      
      if (miss_number==1){
        con=abs(rate-miss_rate*0.01)>0.02
      } else {
        con=((abs(rate-miss_rate*0.01)>0.02) | (sum(miss_each[,1])>sum(miss_each[,2])))
      }
      
    }
    for (ii in 1:miss_number){
      E_cc[dead_id[miss_each[,ii]==0],miss_id[ii]]=NA
    }
    
    miss1=matrix(1,n,1)
    miss1[dead_id]=miss
    miss=miss1
    
    miss_each1=matrix(1,n,miss_number)
    miss_each1[dead_id,]=miss_each
    miss_each=miss_each1
    
    
    
    
  }
  
  print(paste0('Miss rate:',sum(miss==0)/n))
  
  
  ############ KM weight for AFT model ##
  
  if (model=='AFT'){
    Surv.rsp <- Surv(time=Y_cc[,1],event=Y_cc[,2])
    Surv.rsp.new <- Surv(time=Y_test[,1],event=Y_test[,2])
    
    w_full=kmw(Y[,1],Y[,2])
    delta=Y[,2]
    Y=Y[,1]
  } else {
    w_full=NULL
    delta_full=matrix(1,n,1)
  }
  
  
  
  knots=list()
  Boundary.knots=matrix(0,p+q,2)
  for (i in 1:q){
    knots[[i]]=c(0,1)
    Boundary.knots[i,]=c(0,1)
  }
  
  
  para_temp=c(alpha_true,beta_true)
  inter=para_temp==1
  linear_id=data$linear_id==1
  
  
  
  im_time=10
  ########### imputation#########
  fit_impute=boosting_impute(E_cc,im_time,NN,num.knots=c(2),knots=knots,Boundary.knots=Boundary.knots,degree=c(2),'LS',E_type,'nonlinear')$E_impute
  
  
  


  
  
  
  
  
  print('Interaction analysis')
  
  
  fit=list()
  
  para_value=list()
      
  for (i in 1:im_time){
    print(paste('The ',i,'th impute'))
    fit[[i]]=boosting.robust.spline.hier(fit_impute[[i]],G,Y,NN,num.knots=2,knots=knots,Boundary.knots=Boundary.knots,degree=2,'Robust',E_type,'nonlinear',kmweight=w_full)
    para_value[[i]]=variable_return(fit[[i]],p,q)$unique_set
    }
      

  fit_com<-combine_results(para_value,tau,p,q)
  para_hat<-fit_com$para_id
  unique_set<-fit_com$unique_set
  
  y_predict=predict_boosting(E_test,G_test,unique_set,fit[[i]])
  
  
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
    
  estimation_results=estimation_return(unique_set,fit[[i]],p,q)
  
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
  
  
  
