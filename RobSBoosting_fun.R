

library(MASS)
library(pcaPP)
library(splines)
library(Hmisc)


###############################

corr_setting<-function(p,q,corr_str,rho){
  if (corr_str=='AR'){
    sigmaG<-matrix(0,nrow=p,ncol=p)
    diag(sigmaG)<-rep(1,p)
    for(i in 1:p){
      for(j in 1:(i-1)){
        sigmaG[i,j]<-sigmaG[j,i]<-rho^(i-j)
      }
    }
    sigmaE<-matrix(0,nrow=q,ncol=q)
    diag(sigmaE)<-rep(1,q)
    for(i in 1:q){
      for(j in 1:(i-1)){
        sigmaE[i,j]<-sigmaE[j,i]<-rho^(i-j)
      }
    }
  }else{
    sigmaG=matrix(0,p,p)
    diag(sigmaG)=0.5
    for (i in 1:(p-1)){
      sigmaG[i,i+1]=0.33;
    }
    sigmaG=sigmaG+t(sigmaG)
    
    sigmaE=matrix(0,q,q)
    diag(sigmaE)=0.5
    for (i in 1:(q-1)){
      sigmaE[i,i+1]=0.33
    }
    sigmaE=sigmaE+t(sigmaE)
  }
  return(list(sigmaG=sigmaG,sigmaE=sigmaE))
}



simulated_data<-function(n,q,p,G_data_class,model,sigmaG,sigmaE,cc_rate,a=NULL,b=NULL,error_type){
  
  
  G=mvrnorm(n,rep(0,p),sigmaG) #data_class='continuous gene'
  
  
  if (G_data_class=='SNP'){
    for (i in 1:n){
      temp1=quantile(G[i,],0.25)
      temp2=quantile(G[i,],0.75)
      temp3=G[i,]
      G[i,temp3<=temp1]=0;
      G[i,(temp3>temp1)&(temp3<temp2)]=1;
      G[i,temp3>=temp2]=2;
    }
    
    G=scale(G)
    
  }
  
  
  
  
  EC = matrix(runif(n*6,0,1),n,6)
  ED = matrix(rbinom(n*4,1,0.6),n,4)
  
  
 
  EC[,1]=2*sin(EC[,2]*pi*2)+(2*exp(2*EC[,3]-1)-exp(1)+exp(-1))+(-12*EC[,4]*(1+EC[,4])+10)+0.5*ED[,1]+0.5*ED[,2] 
  
  E_coef=matrix(0,q-1,1)
  E_coef[c(1,2,3)]=1
  E_coef[c(6,7)]=0.5
  
  EC[,1]=EC[,1]+rnorm(n,0,1) 
  
  minv=min(EC[,1])
  diffv=max(EC[,1])-min(EC[,1])
  EC[,1]=(EC[,1]-minv)/diffv
  
  
  
  EC_new=EC
  EC_new[,1] = 6*sin(EC[,1]*pi*2)-0.06 #6
  EC_new[,2] = 6*exp(2*EC[,2]-1)-7.05 #6
  alphaC=matrix(0,6,1)
  alphaC[1:2]=1
  
  
  
  EC_new_inter=matrix(0,n,8)
  EC_new_inter[,1]=-2*EC[,1]*(1+EC[,1])+2
  EC_new_inter[,2]=-2*EC[,1]*(1+EC[,1])+2
  EC_new_inter[,3]=-2*EC[,1]*(1+EC[,1])+2
  EC_new_inter[,4] = -2*EC[,1]*(1+EC[,1])+2
  
  EC_new_inter[,5]=-4*(EC[,2])^3+1
  EC_new_inter[,6]=-4*(EC[,2])^3+1
  EC_new_inter[,7]=-4*(EC[,2])^3+1
  EC_new_inter[,8] = -4*(EC[,2])^3+1
  
  a1=1
  a2=1.5
  
  a3=1
  a4=1.5
  
  alphaD=matrix(c(runif(1,a3,a4),0,0,0),4,1)
  
  
  E_new=cbind(EC_new,ED)
  E=cbind(EC,ED)
  a=c(alphaC,alphaD)
  
  
  
  r=matrix(0,p,1)
  r[1:8]=runif(8,a1,a2)
  
  
  beta=matrix(0,q,p)
  
  beta[c(1,2),c(1,2,3,4)]=1
  beta[7,c(5,6,7)]=runif(3,a3,a4)
  
  
  b=matrix(0,q*p+p,1)
  for (j in 1:p) {
    b[(j-1)*(q+1)+1]=r[j]
    b[((j-1)*(q+1)+2):(j*(q+1))]=beta[,j]
  }
  
  
  W=matrix(0,n,p*(q+1))
  W[,seq(from=1,to=p*(q+1),by=(q+1))]=G
  for (i in 1:n){
    temp3=matrix(E[i,],q,1)%*%G[i,]
    temp3[1,1]=EC_new_inter[i,1]*G[i,1]
    temp3[1,2]=EC_new_inter[i,2]*G[i,2]
    temp3[1,3]=EC_new_inter[i,3]*G[i,3]
    temp3[1,4]=EC_new_inter[i,4]*G[i,4]
    temp3[2,1]=EC_new_inter[i,5]*G[i,1]
    temp3[2,2]=EC_new_inter[i,6]*G[i,2]
    temp3[2,3]=EC_new_inter[i,7]*G[i,3]
    temp3[2,4]=EC_new_inter[i,8]*G[i,4]
    
    W[i,setdiff(seq(from=1,to=p*(q+1),by=1),seq(from=1,to=p*(q+1),by=q+1))]=matrix(temp3,p*q,1)
  }
  
  WW=E_new%*%matrix(a,q,1)+W%*%b
  
  E_type=c('EC','EC','EC','EC','EC','EC','ED','ED','ED','ED')
  
  
  
  
  if (error_type==1){
    e=rnorm(n,0,1)
  } else if (error_type==2){
    nn=floor(n*0.9)
    e = as.matrix(c(rnorm(nn), rcauchy(n-nn,0,5)))
  } else if (error_type==3) {
    nn=floor(n*0.9)
    e = as.matrix(c(rnorm(nn), rlnorm(n-nn,0,3)))
  } else if (error_type==4){
    nn=floor(n*0.9)
    temp1=rnorm(n-nn)
    temp2=runif(n-nn,0,1)
    e = as.matrix(c(rnorm(nn), temp1/temp2))
  }
  e = e[sample(n)]
  
  if (model=='linear'){
    YY=WW+e
  } else {
    TT=WW+e
    D=matrix(0,n,1)
    while (abs(1-sum(D)/n-cc_rate*0.01)>0.01){
      if (cc_rate==20) {
        C=runif(n,quantile(TT,0.75),quantile(TT,0.90))
      } else if (cc_rate==40){
        C=runif(n,quantile(TT,0.25),quantile(TT,0.90))
      } 
      Y=pmin(TT,C)
      D=TT<=C+0
    }
    print(paste0('the censoring rate=',1-sum(D)/n))
    YY=matrix(0,n,2)
    YY[,1]=Y
    YY[,2]=D
  }
  
  linear_id=matrix(0,p+q+p*q)
  
  linear_id[7:10]=1
  
  linear_id[(q+1)*(1:p)]=1
  
  linear_id[(q+1)*(1:p)+7]=1
  linear_id[(q+1)*(1:p)+8]=1
  linear_id[(q+1)*(1:p)+9]=1
  linear_id[(q+1)*(1:p)+10]=1
    
    
  para_true=list(alpha=a,beta=b)
 
  
  return(list(G=G,E=E,W=W,para_true=para_true,Y=YY,E_type=E_type,E_new=E_new,E_coef=E_coef,linear_id=linear_id))
  
}

######################

norm2 <- function(a) sqrt(sum(a^2))


######################

basis.matrix<-function(x,num.knots,knots=NULL,Boundary.knots=NULL,degree=3,NorM){
  
  n=length(x)
  
  if (is.null(knots)){
    X = bs(x, df=num.knots+degree+1, intercept=TRUE, degree=degree)
    knots=attr(X,'knots')
    Boundary.knots=attr(X,'Boundary.knots')
  } else {
    X = bs(x, knots=knots, intercept=TRUE, degree=degree,Boundary.knots = Boundary.knots)
    knots=attr(X,'knots')
    Boundary.knots=attr(X,'Boundary.knots')
  }
  
  X=X[,-1]
  X=X-(matrix(1,n,1)%*%NorM)
  
  X=cbind(1,X)
  
  #X = cbind(1,X[,-1])
  
  return(list(X=X,knots=knots,Boundary.knots=Boundary.knots))
}


####################################

design.matrix<-function(x,num.knots,knots,Boundary.knots,degree,v_type,model_type,NorM=NULL){
  ##v_type  EC: continuous E factor, ED: decrete E factor, G: G factor
  ## x EC,ED,G: n*1, G-EC, G-ED: n*2  x[,1]: E factor x[,2]: G factor
  
  
  if (model_type=='linear'){
    
    if ((v_type=='EC') | (v_type=='G') | (v_type=='ED')){
      n=length(x)
      x=matrix(x,n,1)
      X=cbind(1,x)
    } else {
      X=cbind(1,x[,1]*x[,2])
    }
    num.knots=NULL
    knots=NULL
    Boundary.knots=NULL
    
  } else {
    
    if (v_type=='EC'){
      temp=basis.matrix(x,num.knots,knots,Boundary.knots,degree,NorM)
      X=temp$X
      knots=temp$knots
      Boundary.knots=temp$Boundary.knots
    } else if ((v_type=='G') | (v_type=='ED')){
      n=length(x)
      x=matrix(x,n,1)
      X=cbind(1,x)
      num.knots=NULL
      knots=NULL
      Boundary.knots=NULL
    } else if (v_type=='G-EC'){
      temp=basis.matrix(x[,1],num.knots,knots,Boundary.knots,degree,NorM)
      X=temp$X
      d_temp=dim(X)[2]
      X[,2:d_temp]=X[,2:d_temp]*(x[,2]%*%matrix(1,1,d_temp-1))
      knots=temp$knots
      Boundary.knots=temp$Boundary.knots
    } else if (v_type=='G-ED'){
      n=dim(x)[1]
      X=cbind(1,x[,1]*x[,2])
      num.knots=NULL
      knots=NULL
      Boundary.knots=NULL
    }
  }
  return(list(X=X,knots=knots,Boundary.knots=Boundary.knots,degree=degree))
}


##########################

pen.ls <- function(y, X,kmweight=NULL)
{
  
  if (is.null(kmweight)) {
    beta.ls=matrix(0,dim(X)[2],1)
    temp=solve( t(X[,colSums(X)!=0]) %*% X[,colSums(X)!=0])
    beta.ls[colSums(X)!=0] <- as.vector(temp %*% t(X[,colSums(X)!=0]) %*% y )
  } else {
    beta.ls=matrix(0,dim(X)[2],1)
    Xtemp=X*kmweight
    ctemp= solve(t(Xtemp[,colSums(Xtemp)!=0])%*%X[,colSums(Xtemp)!=0])%*%(t(Xtemp[,colSums(Xtemp)!=0]))  
    beta.ls[colSums(Xtemp)!=0] <-ctemp%*% y
  }
  y_hat=X%*%beta.ls
  
  
  return(list(beta=beta.ls,y_hat=y_hat))
}




pen.m<- function(y,X,r,kmweight=NULL) 
{
  n=length(y)
  rr=psi.huber(r,k=1.345*mad(r))
  W=rr
  if (!is.null(kmweight)){
    W=W*kmweight
  }
  
  estimates=matrix(0,dim(X)[2],1)
  
  Xtemp=X*W
  ctemp= solve(t(Xtemp[,colSums(Xtemp)!=0])%*%X[,colSums(Xtemp)!=0])%*%(t(Xtemp[,colSums(Xtemp)!=0]))
  mbeta <-ctemp%*% y
  y_hat=X[,colSums(Xtemp)!=0]%*%mbeta
  estimates[colSums(Xtemp)!=0]=mbeta
  
  return(list(y_hat=y_hat,estimates=estimates,scale=scale))
 
}

####################################


robust.spline<-function(designmat,y,r,kmweight=NULL){
  
  
  
  X=designmat$X
  knots=designmat$knots
  Boundary.knots=designmat$Boundary.knots
  degree=designmat$degree
  
  fit<-pen.m(y, X,r,kmweight) 
  return(list(estimates=fit$estimates,y_hat=fit$y_hat,knots=knots,Boundary.knots=Boundary.knots,degree=degree,X=X))
}


#######################

ls.spline<-function(designmat,y,kmweight=NULL){
  
  
  X=designmat$X
  knots=designmat$knots
  Boundary.knots=designmat$Boundary.knots
  degree=designmat$degree
  
  fit<-pen.ls(y, X,kmweight) 
  
  return(list(estimates=fit$beta,y_hat=fit$y_hat,knots=knots,Boundary.knots=Boundary.knots,degree=degree,X=X))
}


#########################

boosting.robust.spline<-function(x,y,loop_time,num.knots=NULL,knots=NULL,Boundary.knots=NULL,degree=1,Method,E_type,model_type,kmweight=NULL,v=0.1){
  
  
  if (is.matrix(x)){
    n=dim(x)[1]
    p=dim(x)[2]
  } else {
    p=length(x)
    n=length(x[[1]])
  }
  
 
  if (!is.null(knots)){
    ss<-seq(from=1/(num.knots+1),to=num.knots/(num.knots+1),by=1/(num.knots+1))
    knots_temp=vector('list',p)
    for (ii in 1:q){
      knots_temp[[ii]]=quantile(knots[[ii]],ss)
    }
    knots=knots_temp
    xx_dots=seq(from=0,to=1,by=0.0001)
    NorM = bs(xx_dots, knots=knots[[1]],degree=degree,Boundary.knots = Boundary.knots[1,])
    NorM=colMeans(NorM)
  } else {
    xx_dots=seq(from=0,to=1,by=0.0001)
    NorM = bs(xx_dots, df=num.knots+degree,degree=degree)
    NorM=colMeans(NorM)
  }
  
  
 
  u=y
  f=0
  
  result=list()
  variable=matrix(0,loop_time,1)
  BIC=matrix(-1e+20,loop_time,1)
  loglike=matrix(0,loop_time,1)
  t=1
  f_temp=f
  
  
  designmat=vector('list',p)
  v_type=E_type
  
  
  
  if (is.matrix(x)){
    
    if (is.null(knots)){
      
      for (j in 1:p){
        designmat[[j]]=design.matrix(x[,j],num.knots,knots,Boundary.knots,degree,v_type[j],model_type,NorM)
      }
    } else {
      for (j in 1:p){
        designmat[[j]]=design.matrix(x[,j],num.knots,knots[[j]],Boundary.knots[j,],degree,v_type[j],model_type,NorM)
      }
    }
  } else {
    if (is.null(knots)){
      for (j in 1:p){
        designmat[[j]]=design.matrix(x[[j]],num.knots,knots,Boundary.knots,degree,v_type[j],model_type,NorM)
      }
    } else {
      for (j in 1:p){
        designmat[[j]]=design.matrix(x[[j]],num.knots,knots[[j]],Boundary.knots[j,],degree,v_type[j],model_type,NorM)
      }
    }
    
    
    
    
  }
  
  
  hatmat_set=list()
  
  
  
  for (t in 1:loop_time){
    res=matrix(1e+20,p,1)
    for (j in 1:p){
      if (Method=='Robust'){
        temp=robust.spline(designmat[[j]],u,y-f,kmweight)
      } else if (Method=='LS') {
        temp=ls.spline(designmat[[j]],u,kmweight)
      }
      
      y_predict=temp$y_hat
      f1=y_predict+f_temp
      hatmat_set[[j]]=temp
      
      
      if (Method=='Robust'){
        ptemp=dim(temp$X)[2]
        if (ptemp==2){
          sdd1=sqrt(mean((temp$X[,2:ptemp]-mean(temp$X[,2:ptemp]))^2))
          sdd2=qn(u)
          RSS=n*(sdd2^2-(temp$estimate[2]*sdd1)^2)
        } else {
          xtemp=temp$X[,2:ptemp]
          xxtemp=xtemp-matrix(1,n,1)%*%colMeans(xtemp)
          robustcov=(t(xxtemp)%*%xxtemp)/n
          beta_temp=matrix(temp$estimates[2:ptemp],1,ptemp-1)
          RSS=n*(qn(u)^2-beta_temp%*%robustcov%*%t(beta_temp))
        }
      } else if (Method=='LS') {
        RSS=sum((y-f1)^2)
      }
      if (RSS<0){
        RSS=1e+100
      }
      df=length(unique(c(variable[1:(t-1)],j)))
      res[j]=log(RSS)+df*log(max(n,p))/n
      
    }
    
    id=which.min(res)
    temp=hatmat_set[[id]]
    y_predict=temp$y_hat
    coef=temp$estimates
    result[[t]]=temp
    variable[t]=id
    f=f+v*y_predict
    f_temp=f
    u=u-v*y_predict
    
   
    
    
    if (Method=='Robust') {
      absy=abs(y-f)
      MAD_v=mad(y-f)
      
      c = 1.345*MAD_v
      Huber=absy^2
      Huber[absy>c]=2*c*(absy[absy>c]-c/2)
      RSS=sum(Huber)
    } else if (Method=='LS') {
      RSS=sum((y-f)^2)
    }
    
    if (RSS==0){
      RSS=1e-200
    }
    df=length(unique(variable[1:t]))
    BIC[t]=log(RSS) + df*log(n)/n
    if (t>1000) {
      if (abs((BIC[t]-BIC[t-1])/BIC[t-1])<1e-8){
        break
      }
    }
    t=t+1
  }
  
  BIC=BIC[1:t]
  variable=variable[1:t]
  loglike=loglike[1:t]
  
  id=which.min(BIC)
  output=list(spline_result=result,BIC=BIC,variable=variable,id=id,max_t=t,model_type=model_type,degree=degree,v=v,NorM=NorM)
  
  return(output)
  
  
}


###############################


boosting.robust.spline.hier<-function(E,G,y,loop_time,num.knots=NULL,knots=NULL,Boundary.knots=NULL,degree=1,Method,E_type,model_type,kmweight=NULL,v=0.1){
  
  
  
  n=dim(E)[1]
  p=dim(G)[2]
  q=dim(E)[2]
  
  if (!is.null(knots)){
    ss<-seq(from=1/(num.knots+1),to=num.knots/(num.knots+1),by=1/(num.knots+1))
    knots_temp=vector('list',p+q)
    for (ii in 1:q){
      knots_temp[[ii]]=quantile(knots[[ii]],ss)
    }
    knots=knots_temp
    xx_dots=seq(from=0,to=1,by=0.0001)
    NorM = bs(xx_dots, knots=knots[[1]],degree=degree,Boundary.knots = Boundary.knots[1,])
    NorM=colMeans(NorM)
  } else {
    xx_dots=seq(from=0,to=1,by=0.0001)
    NorM = bs(xx_dots, df=num.knots+degree,degree=degree)
    NorM=colMeans(NorM)
  }
  
  
  
  variable<-BIC<-loglike<-c()
  result=list()
  
  variable_pair=matrix(0,p+q,2)
  variable_pair[,1]=c(1:q,rep(0,p))
  variable_pair[,2]=c(rep(0,q),1:p)
  
  colnames(variable_pair)<-c('E','G')
  
  GE=cbind(E,G)
  
  v_type=c(E_type,rep('G',p))
  
  designmat=vector('list',dim(GE)[2])
  
  
  if (is.null(knots)){
    
    for (j in 1:dim(GE)[2]){
      designmat[[j]]=design.matrix(GE[,j],num.knots,knots,Boundary.knots,degree,v_type[j],model_type,NorM)
    }
  } else {
    for (j in 1:dim(GE)[2]){
      designmat[[j]]=design.matrix(GE[,j],num.knots,knots[[j]],Boundary.knots[j,],degree,v_type[j],model_type,NorM)
    }
  }
  
 
  f=0
  u=y-f
  
  
  t=1
  f_temp=f

  
  while (t<=loop_time) {
    pp=length(designmat)
    res=matrix(1e+20,pp,1)
    hatmat_set=vector('list',pp)
    
  
    for (j in 1:pp){
      
      if (Method=='Robust'){
        temp=robust.spline(designmat[[j]],u,y-f,kmweight)
      } else if (Method=='LS') {
        temp=ls.spline(designmat[[j]],u,kmweight)
      }
     
      y_predict=temp$y_hat
      f1=y_predict+f_temp
      hatmat_set[[j]]=temp
      
      if (Method=='Robust'){
        ptemp=dim(temp$X)[2]
        
        if (ptemp==2){
          sdd1=sqrt(mean((temp$X[,2:ptemp]-mean(temp$X[,2:ptemp]))^2))
          sdd2=qn(u)
          
          RSS=n*(sdd2^2-(temp$estimate[2]*sdd1)^2)
          
        } else {
          xtemp=temp$X[,2:ptemp]
          xxtemp=xtemp-matrix(1,n,1)%*%colMeans(xtemp)
          robustcov=(t(xxtemp)%*%xxtemp)/n
          beta_temp=matrix(temp$estimates[2:ptemp],1,ptemp-1)
          RSS=n*(qn(u)^2-beta_temp%*%robustcov%*%t(beta_temp))
        }
      } else if (Method=='LS') {
        RSS=sum((y-f1)^2)
      }
      if (RSS<0){
        RSS=1e+100
      }
      
      df=length(unique(c(variable[1:(t-1)],j)))
      
      res[j]=log(RSS)+df*log(n)/n
      
    }
    
    
    
    id=which.min(res)
    temp=hatmat_set[[id]]
    y_predict=temp$y_hat
    coef=temp$coef
    
    result[[t]]=temp
    variable[t]=id
    f=f+v*y_predict
    f_temp=f
    
    
    u=y-f
    
    id_unique=unique(variable[1:t])
    
    
    if (length(id_unique)>1){ #### have at least two main effects
      if (variable[t]<=p+q){ #### the selected variable is a main effect
        if (!is.element(variable[t],variable[1:(t-1)])){ # the selected variable is a new main effect
          id_E=id_unique[which(id_unique<=q)] 
          id_G=id_unique[which((id_unique<=(p+q)) & (id_unique>q))]-q
          if ((length(id_E)>0) & (length(id_G)>0)) {# have at least one E and one G effects
            if (variable[t]<=q) {
              
              temp_pair=matrix(0,length(id_G),2)
              temp_pair[,1]=rep(variable[t],length(id_G))
              temp_pair[,2]=id_G
              vtype_add=matrix(0,1,length(id_G))
              ppp=pp+length(id_G)
              for (ii in 1:length(id_G)){
                if (E_type[variable[t]]=='ED'){
                  vtype_add[ii]='G-ED'
                } else {
                  vtype_add[ii]='G-EC'
                }
                temp_data=cbind(E[,variable[t]],G[,id_G[ii]])
                
                
                
                if (is.null(knots)){
                  designmat[[pp+ii]]=design.matrix(temp_data,num.knots,knots,Boundary.knots,degree,vtype_add[ii],model_type,NorM)
                } else {
                  designmat[[pp+ii]]=design.matrix(temp_data,num.knots,knots[[variable[t]]],Boundary.knots[variable[t],],degree,vtype_add[ii],model_type,NorM) 
                }
                
              }
            } else {
              
              temp_pair=matrix(0,length(id_E),2)
              temp_pair[,1]=id_E
              temp_pair[,2]=rep(variable[t]-q,length(id_E))
              vtype_add=matrix(0,1,length(id_E))
              ppp=pp+length(id_E)
              for (ii in 1:length(id_E)){
                if (E_type[id_E[ii]]=='ED'){
                  vtype_add[ii]='G-ED'
                } else {
                  vtype_add[ii]='G-EC'
                }
                temp_data=cbind(E[,id_E[ii]],G[,variable[t]-q])
                if (is.null(knots)){
                  designmat[[pp+ii]]=design.matrix(temp_data,num.knots,knots,Boundary.knots,degree,vtype_add[ii],model_type,NorM)
                } else {
                  
                  designmat[[pp+ii]]=design.matrix(temp_data,num.knots,knots[[id_E[ii]]],Boundary.knots[id_E[ii],],degree,vtype_add[ii],model_type,NorM) 
                  
                  
                }
              }
            }
            
            v_type=c(v_type,vtype_add)
            
            variable_pair=rbind(variable_pair,temp_pair)
          }
        }
      }
    }
    
   
    
    if (Method=='Robust') {
      absy=abs(y-f)
      MAD_v=mad(y-f)
      c = 1.345*MAD_v
      Huber=absy^2
      Huber[absy>c]=2*c*(absy[absy>c]-c/2)
      RSS=sum(Huber)
    } else if (Method=='LS') {
      RSS=sum((y-f)^2)
    }
    
    df=length(unique(variable[1:t]))
    BIC[t]=log(RSS)+df*log(n)/(n)
    
    
    if (t>1000) {
      if (abs((BIC[t]-BIC[t-1])/BIC[t-1])<1e-8){
        break
      }
    }
    t=t+1
  }
  
  
  
  
  
  
  
  
  
  BIC=BIC[1:t]
  variable=variable[1:t]
  loglike=loglike[1:t]
  
  id=which.min(BIC)
  
  output=list(spline_result=result,BIC=BIC,variable=variable,id=id,GE=GE,variable_pair=variable_pair,max_t=t,v_type=v_type,model_type=model_type,degree=degree,v=v,NorM=NorM)
  return(output)
}


variable_return<-function(fit,p,q){
  
  variable_pair=unique(fit$variable_pair[fit$variable[1:fit$id],],MARGIN=1)
  if (!is.matrix(variable_pair)){
    variable_pair=matrix(variable_pair,1,2)
  }
  
  alpha0=matrix(0,q,1)
  G_main=matrix(0,p,1)
  GE_interaction=matrix(0,q,p)
  temp=variable_pair[which(variable_pair[,2]==0),1]
  alpha0[temp]=1
  temp=variable_pair[which(variable_pair[,1]==0),2]
  G_main[temp]=1
  
  temp=variable_pair[((variable_pair[,1]!=0) & (variable_pair[,2]!=0)),]
  
  GE_interaction[temp]=1
  
  beta0=matrix(0,p+p*q,1)
  
  for (j in 1:p) {
    beta0[(j-1)*(q+1)+1]=G_main[j]
    beta0[((j-1)*(q+1)+2):(j*(q+1))]=GE_interaction[,j]
  }
  
  
  
  
  unique_temp=unique(fit$variable[1:fit$id])
  unique_variable=fit$variable_pair[unique_temp,]
  unique_coef=vector('list',length(unique_temp))
  unique_knots=vector('list',length(unique_temp))
  unique_Boundary.knots=vector('list',length(unique_temp))
  for (i in 1:length(unique_temp)){
    unique_coef[[i]]=0
  }
  
  
  v=fit$v
  a=0
  for (i in 1:fit$id){
    id_temp=which(unique_temp==fit$variable[i])
    unique_coef[[id_temp]]=unique_coef[[id_temp]]+fit$spline_result[[i]]$estimates[-1]*v
    a=a+fit$spline_result[[i]]$estimates[1]*v
    if ((!is.null(fit$spline_result[[i]]$knots))){
      unique_knots[[id_temp]]=fit$spline_result[[i]]$knots
      unique_Boundary.knots[[id_temp]]=fit$spline_result[[i]]$Boundary.knots
    }
  }
  
  
  
  unique_vtype=fit$v_type[unique_temp] 
  
  
  
  unique_set=list(intercept=a,unique_variable=unique_variable,unique_coef=unique_coef,unique_knots=unique_knots,unique_Boundary.knots=unique_Boundary.knots,unique_vtype=unique_vtype)
  
  
  return(list(para_id=c(alpha0,beta0),unique_set=unique_set))
}


estimation_return<-function(unique_set,fit,p,q){
  id=length(unique_set$unique_coef)
  unique_variable=unique_set$unique_variable
  unique_coef=unique_set$unique_coef
  unique_knots=unique_set$unique_knots
  unique_Boundary.knots=unique_set$unique_Boundary.knots
  unique_vtype=unique_set$unique_vtype
  
  if (id==1){
    unique_variable=matrix(unique_variable,1,2)
  }
  
  
  xx_dot=seq(from=0,to=1,by=0.0001)
  
  estimation_results=vector('list',p+q+p*q) ## E+G+E*G
  
  
  
  
  id_temp1=which(unique_set$unique_vtype=='EC')
  id_temp2=unique_set$unique_variable[id_temp1,1]+unique_set$unique_variable[id_temp1,2]*(q+1)
  
  for (i in 1:length(id_temp1)){
    if (fit$model_type=='nonlinear') {
      xx=bs(xx_dot, knots=unique_set$unique_knots[[id_temp1[i]]], intercept=TRUE, degree=fit$degree,Boundary.knots = unique_set$unique_Boundary.knots[[id_temp1[i]]])
      xx=xx[,-1]
      xx=xx-(matrix(1,length(xx_dot),1)%*%fit$NorM)
      bs_predict=xx%*%unique_set$unique_coef[[id_temp1[i]]]
    } else {
      bs_predict=matrix(unique_set$unique_coef[[id_temp1[i]]],length(xx_dot),1) 
    }
    estimation_results[[id_temp2[i]]]=bs_predict
  }
  
  id_temp1=which(unique_set$unique_vtype=='G-EC')
  id_temp2=unique_set$unique_variable[id_temp1,1]+unique_set$unique_variable[id_temp1,2]*(q+1)
  id_temp3=unique_set$unique_variable[id_temp1,2]
  
  for (i in 1:length(id_temp1)){
    if (fit$model_type=='nonlinear') {
      xx=bs(xx_dot, knots=unique_set$unique_knots[[id_temp1[i]]], intercept=TRUE, degree=fit$degree,Boundary.knots = unique_set$unique_Boundary.knots[[id_temp1[i]]])
      xx=xx[,-1]
      xx=xx-(matrix(1,length(xx_dot),1)%*%fit$NorM)
      bs_predict=xx%*%unique_set$unique_coef[[id_temp1[i]]]
    } else {
      bs_predict=matrix(unique_set$unique_coef[[id_temp1[i]]],length(xx_dot),1) 
    }
    estimation_results[[id_temp2[i]]]=bs_predict
  }
  
  
  id_temp1=which(((unique_set$unique_vtype=='ED') | (unique_set$unique_vtype=='G') | (unique_set$unique_vtype=='G-ED'))!=0)
  id_temp2=unique_set$unique_variable[id_temp1,1]+unique_set$unique_variable[id_temp1,2]*(q+1)
  
  
  for (i in 1:length(id_temp1)){
    estimation_results[[id_temp2[i]]]=unique_set$unique_coef[[id_temp1[i]]]
  }
  
  intercept=unique_set$intercept
  return(list(intercept=intercept,estimation_results=estimation_results))
}



predict_boosting<-function(E_test,G_test,unique_set,fit){
  
  id=length(unique_set$unique_coef)
  unique_variable=unique_set$unique_variable
  unique_coef=unique_set$unique_coef
  unique_knots=unique_set$unique_knots
  unique_Boundary.knots=unique_set$unique_Boundary.knots
  unique_vtype=unique_set$unique_vtype
  
  if (id==1){
    unique_variable=matrix(unique_variable,1,2)
  }
  
  y_predict=0
  
  for (t in 1:id){
    if (unique_variable[t,1]==0){
      x_temp=G_test[,unique_variable[t,2]]
    } else if (unique_variable[t,2]==0){
      x_temp=E_test[,unique_variable[t,1]]
    } else {
      x_temp=cbind(E_test[,unique_variable[t,1]],G_test[,unique_variable[t,2]])
    }
    
    N=design.matrix(x_temp,knots=unique_knots[[t]],Boundary.knots=unique_Boundary.knots[[t]],degree=fit$degree,v_type=unique_vtype[t],model_type=fit$model_type,NorM=fit$NorM)$X
    temp=as.matrix(N[,-1])%*%unique_coef[[t]]
    y_predict=y_predict+temp
  }
  y_predict=y_predict+unique_set$intercept
  return(y_predict)
}


predict_boosting2<-function(fit,x_new,v_type){
  
  y_predict=0
  id=fit$id
  v=fit$v
  for (t in 1:id){
    knots=fit$spline_result[[t]]$knots
    coef=fit$spline_result[[t]]$estimates
    degree=fit$spline_result[[t]]$degree
    Boundary.knots=fit$spline_result[[t]]$Boundary.knots
    if (is.list(x_new)){
      x_temp=x_new[[fit$variable[t]]]
    } else {
      x_temp=x_new[,fit$variable[t]]
    }
    N=design.matrix(x_temp,knots=knots,Boundary.knots=Boundary.knots,degree=degree,v_type=v_type[fit$variable[t]],model_type=fit$model_type,NorM=fit$NorM)$X
    temp=N%*%coef
    y_predict=y_predict+v*temp
  }
  return(y_predict)
}


boosting_impute<-function(E_cc,im_time,loop_time,num.knots=NULL,knots=NULL,Boundary.knots=NULL,degree=1,Method,E_type,model_type){
  
  
  
  n=dim(E_cc)[1]
  q=dim(E_cc)[2]
  
  
  miss_id=which(colSums(is.na(E_cc))!=0)
  
  p1=1
  p2=q-1
  
  
  
  E_impute=list()
  
  fit_set=list()
  
  if (length(miss_id)==1){
    
    miss=matrix(0,n,1)
    miss[rowSums(is.na(E_cc))==0]=1
    
    nonmiss_id=which(colSums(is.na(E_cc))==0)
    n_m=sum(miss==0)
    n_c=sum(miss==1)
    
    E_miss_obs=matrix(E_cc[miss==1,miss_id],n_c,p1) # observed samples for the miss variable  y_train for boosting spline
    E_miss_miss=matrix(E_cc[miss==0,miss_id],n_m,p1) # missing samples for the miss variable  y_test for boosting spline
    
    E_obs_obs=matrix(E_cc[miss==1,nonmiss_id],n_c,p2) # observed samples for the nonmiss variable x_train for boosting spline
    E_obs_miss=matrix(E_cc[miss==0,nonmiss_id],n_m,p2) # missing samples for the nonmiss variable  x_test for boosting spline
    
    fit<-boosting.robust.spline(E_obs_obs,E_miss_obs,loop_time,num.knots,knots,Boundary.knots,degree,Method,E_type[nonmiss_id],model_type)
    
    fit_set[[1]]=fit
    
    im_mean=predict_boosting2(fit,E_obs_miss,E_type[nonmiss_id])
    
    
    nnn=n_c-ceiling(n_c/2)
    
    err=matrix(0,nnn,1)
    
    repeat_time=10
    
    
    for (rr in 1:repeat_time){
      D1=sample(n_c)[1:ceiling(n_c/2)]
      D2=setdiff(1:n_c,D1)
      fit1<-boosting.robust.spline(E_obs_obs[D1,],E_miss_obs[D1,],loop_time,num.knots,knots,Boundary.knots,degree,Method,E_type[nonmiss_id],model_type)
      predict_cv=predict_boosting2(fit1,E_obs_obs[D2,],E_type[nonmiss_id])
      err[((rr-1)*nnn+1):(rr*nnn)]=E_miss_obs[D2,]-predict_cv
    }
    
    im_value=matrix(0,n_m,im_time)
    
    
    
    for (ii in 1:im_time){
      idd=1:n_m
      temp=matrix(0,1,n_m)
      while (length(idd)>0){
        noise_random=emprand(err,length(idd))
        temp[idd]=im_mean[idd]+noise_random
        idd=which((temp>1) | (temp<0))
      }
      im_value[,ii]=temp
    }
    
    for (i in 1:im_time){
      temp=E_cc
      temp[miss==0,miss_id]=im_value[,i]
      E_impute[[i]]=temp
    }
  } else {
    
    
    for (i in 1:im_time){
      E_impute[[i]]=E_cc
    }
    
    miss_id_s=miss_id
    
    miss_each=matrix(0,n,q)
    miss_each[!is.na(E_cc)]=1
    
    miss_rate_each=colSums(miss_each[,miss_id])
    miss_rank=sort(miss_rate_each,index.return=T,decreasing=T)$ix
    
    
    
    for (mm in 1:im_time){
      
      for (ii in 1:length(miss_id_s)) {
        
        ii=miss_rank[ii]
        
        E_temp=E_impute[[mm]]
        
        miss=matrix(0,n,1)
        miss[rowSums(is.na(E_temp))==0]=1
        
        nonmiss_id=which(colSums(is.na(E_temp))==0)
        nn_c=sum(miss==1)
        
        
        nn_m=sum(miss_each[,ii]==0)
        
        p2=length(nonmiss_id)
        
        E_obs_obs=matrix(E_temp[miss==1,nonmiss_id],nn_c,p2) # observed samples for the nonmiss variable x_train for boosting spline
        
        E_obs_miss=matrix(E_temp[miss_each[,ii]==0,nonmiss_id],nn_m,p2) # missing samples for the nonmiss variable  x_test for boosting spline
        E_miss_miss=matrix(E_temp[miss_each[,ii]==0,ii],nn_m,1) # missing samples for the miss variable  y_test for boosting spline
        E_miss_obs=matrix(E_temp[miss==1,ii],nn_c,p1) # observed samples for the miss variable  y_train for boosting spline
        
        
        fit<-boosting.robust.spline(E_obs_obs,E_miss_obs,loop_time,num.knots,knots,Boundary.knots,degree,Method,E_type[nonmiss_id],model_type)
        
        fit_set[[ii]]=fit
        im_mean=predict_boosting2(fit,E_obs_miss,E_type[nonmiss_id])
        
        
        nnn=nn_c-ceiling(nn_c/2)
        
        err=matrix(0,nnn,1)
        
        repeat_time=10
        
        
        for (rr in 1:repeat_time){
          D1=sample(nn_c)[1:ceiling(nn_c/2)]
          D2=setdiff(1:nn_c,D1)
          fit1<-boosting.robust.spline(E_obs_obs[D1,],E_miss_obs[D1,],loop_time,num.knots,knots,Boundary.knots,degree,Method,E_type[nonmiss_id],model_type)
          predict_cv=predict_boosting2(fit1,E_obs_obs[D2,],E_type[nonmiss_id])
          err[((rr-1)*nnn+1):(rr*nnn)]=E_miss_obs[D2,]-predict_cv
        }
        
        
        idd=1:nn_m
        temp=matrix(0,1,nn_m)
        while (length(idd)>0){
          noise_random=emprand(err,length(idd))
          temp[idd]=im_mean[idd]+noise_random
          idd=which((temp>1) | (temp<0))
        }
        im_value=temp
        
        
        
        if (E_type[ii]=='ED'){
          im_value[im_value>=0.5]=1
          im_value[im_value<0.5]=0
        }
        
        temp=E_impute[[mm]]
        temp[miss_each[,ii]==0,ii]=im_value
        E_impute[[mm]]=temp
      }
    }
  }
  return(list(E_impute=E_impute,fit_set=fit_set))
}




combine_results<-function(para_value,tau,p,q){
  
  
  
  im_time=length(para_value)
  
  for (i in 1:im_time){
    if (i==1){
      unique_variable=para_value[[i]]$unique_variable
      unique_knots=para_value[[i]]$unique_knots
      unique_Boundary.knots=para_value[[i]]$unique_Boundary.knots
      intercept=para_value[[i]]$intercept
      unique_vtype=para_value[[i]]$unique_vtype
    } else {
      unique_variable=rbind(unique_variable,para_value[[i]]$unique_variable)
      unique_knots=c(unique_knots,para_value[[i]]$unique_knots)
      unique_Boundary.knots=c(unique_Boundary.knots,para_value[[i]]$unique_Boundary.knots)
      unique_vtype=c(unique_vtype,para_value[[i]]$unique_vtype)
      intercept=intercept+para_value[[i]]$intercept
    }
  }
  
  intercept=intercept/im_time
  
  dp=duplicated(unique_variable,MARGIN=1)
  unique_variable=unique_variable[!dp,]
  unique_knots=unique_knots[!dp]
  unique_Boundary.knots=unique_Boundary.knots[!dp]
  unique_vtype=unique_vtype[!dp]
  unique_coef_set=list()
  
  if (!is.matrix(unique_variable)){
    unique_variable=matrix(unique_variable,1,2)
  }
  
  
  
  for (i in 1:dim(unique_variable)[1]){
    unique_coef_set[[i]]=0
  }
  unique_coef_id=matrix(0,sum(!dp),im_time)
  for (ii in 1:dim(unique_variable)[1]){
    for (i in 1:im_time){
      temp=which.column(unique_variable[ii,],para_value[[i]]$unique_variable)
      if (length(temp)!=0){
        unique_coef_set[[ii]]=unique_coef_set[[ii]]+para_value[[i]]$unique_coef[[temp]]/im_time
        unique_coef_id[ii,i]=1
      }
    }
  }
  
  id=which(rowMeans(unique_coef_id)>tau)
  
  unique_variable=unique_variable[id,]
  unique_knots=unique_knots[id]
  unique_Boundary.knots=unique_Boundary.knots[id]
  unique_vtype=unique_vtype[id]
  unique_coef=unique_coef_set[id]
  
  if (!is.matrix(unique_variable)){
    unique_variable=matrix(unique_variable,1,2)
  }
  
  
  variable_pair=unique_variable
  if (!is.matrix(variable_pair)){
    variable_pair=matrix(variable_pair,1,2)
  }
  
  alpha0=matrix(0,q,1)
  G_main=matrix(0,p,1)
  GE_interaction=matrix(0,q,p)
  temp=variable_pair[which(variable_pair[,2]==0),1]
  alpha0[temp]=1
  temp=variable_pair[which(variable_pair[,1]==0),2]
  G_main[temp]=1
  
  temp=variable_pair[((variable_pair[,1]!=0) & (variable_pair[,2]!=0)),]
  
  GE_interaction[temp]=1
  
  beta0=matrix(0,p+p*q,1)
  
  for (j in 1:p) {
    beta0[(j-1)*(q+1)+1]=G_main[j]
    beta0[((j-1)*(q+1)+2):(j*(q+1))]=GE_interaction[,j]
  }
  
  unique_set=list(intercept=intercept,unique_variable=unique_variable,unique_coef=unique_coef,unique_knots=unique_knots,unique_Boundary.knots=unique_Boundary.knots,unique_vtype=unique_vtype)
  
  
  return(list(para_id=c(alpha0,beta0),unique_set=unique_set))
  
}





which.column<-function(a,b){
  
  if (!is.matrix(b)){
    b=matrix(b,1,2)
  }
  
  id1=is.element(b[,1],a[1])
  id2=is.element(b[,2],a[2])
  id=which(id1*id2==1)
  
  return(id)
}








GetFPTP<-function(theta,theta_hat){
  # to get TNR (True Negative Rate ) and TPR (True Positive Rate) 
  thea = abs(theta) > 0   # transform coefficients to binary values
  thea_hat = abs(theta_hat) > 1e-8  # convert estimated coefficients to binary values
  A = sum((!thea)*(!thea_hat))  # A: TN
  B = sum((!thea)*thea_hat)   # B: FP
  C = sum(thea*(!thea_hat))   # C: FN
  D = sum(thea*thea_hat)    # D: TP
  TPR = D/(D+C)    # TPR=TP/(TP+FN)  true positive rate (TPR) sensitivity
  FPR = A/(B+A)    # TNR=TN/(TN+FP)  true negative rate    specificity 
  result=list(TPR= TPR, FPR = FPR,TP=D,FP=B)
  return(result)
}


emprand<-function(dist,n){
  
  
  x = c(dist)
  #Remove missing observations indicated by NaN's.
  id = !is.na(x)
  x = x[id]
  
  # Compute empirical cumulative distribution function (cdf)
  xlen = length(x)
  x = sort(x)
  p = seq(from=1,to=xlen,by=1)
  p = p/xlen   
  
  # Generate uniform random number between 0 and 1
  ur =  runif(n,0,1)
  
  # Interpolate ur from empirical cdf and extraplolate for out of range  values.
  xr = approxExtrap(p,x,ur,rule=2)$y
  return(xr)
}  







kmw<-function(y,delta){
  y_s=y
  delta_s=delta
  kmweight<-c()
  nw<-length(y)
  
  comb<-cbind(y,delta)
  oo=order(y)
  ocomb<-comb[oo,]
  y<-ocomb[,1]
  delta<-ocomb[,2]
  kmweight[1]<-delta[1]/nw
  for(ind in 2:nw){
    tmp<-c()
    for(ind2 in 1:(ind-1)){
      tmp[ind2]<-((nw-ind2)/(nw-ind2+1))^delta[ind2]
    }
    kmweight[ind]<-delta[ind]/(nw-ind+1)*prod(tmp)
  }
  kmweight1=matrix(0,nw,0)
  kmweight1[oo]=kmweight
  return(kmweight=kmweight1)
}



