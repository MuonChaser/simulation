# source("CADEreg_new.R")
source("assign.R")
library(experiment)
library(ggplot2)

J=10
n.avg=20
source("gen_data.R")

r <- function(data){
if(!is.factor(data$id)){stop('The cluster_id should be a factor variable.')}
  cluster.id=unique(data$id)	
  n.cluster=length(cluster.id)	
  Z=vector("list", n.cluster) 	
  D=vector("list", n.cluster) 
  Y=vector("list", n.cluster) 
  B=vector("list", n.cluster)
  A=rep(0,n.cluster)
  for (i in 1:n.cluster){
    Z[[i]]=as.numeric(data$Z[data$id==cluster.id[i]])
    D[[i]]=as.numeric(data$D[data$id==cluster.id[i]])
    Y[[i]]=data$Y[data$id==cluster.id[i]]
    B[[i]]=as.numeric(data$B[data$id==cluster.id[i]])
    if (length(unique(data$A[data$id==cluster.id[i]]))!=1){
      stop( paste0('The assignment mechanism in cluster ',i,' should be the same.'))
    }
    A[i]=data$A[data$id==cluster.id[i]][1]
  }
  
	
  n=sapply(Z,length)
  J=length(n)
  # weights 
  W=sum(n)
  J1=sum(A)
  J0=J-J1
  n1=sapply(Z,sum)
  n0=n-n1
  # K: number of strata
  K=max(unlist(B))
  # nk[k]: total number of individuals in stratum k
  nk=rep(0,K)
  # njk[j,k]: number of individuals in stratum k within cluster j
  njk=matrix(0, nrow=J, ncol=K)
  for (j in 1:J){
    for (kk in 1:K){
      njk[j,kk]=sum(B[[j]]==kk)
      nk[kk]=nk[kk]+njk[j,kk]
    }
  }
  N.total=sum(n)

  index.l=rep(1,J)
  index.r=rep(1,J)   
  for(j in 2:J){
    index.l[j]=1+sum(n[1:(j-1)])
  } 
  for(j in 1:J){
    index.r[j]=sum(n[1:j])
  }     
  for (j in 1:J){
    index=index.l[j]:index.r[j]
    W[index]=ifelse(A[j]==1, 1/J1,1/(J-J1)) *   ifelse(Z[[j]]==1, 1/n1[j],1/(n0[j]))
  }
  ## Design matrix in the first stage
  A.reg= rep(0,sum(n))
  Z.reg=rep(0,sum(n))
  D.reg=rep(0,sum(n))
  Y.reg=rep(0,sum(n))
  Pi.full=matrix(0, nrow=sum(n), ncol=K)
  for (j in 1:J){
    index=index.l[j]:index.r[j]
    A.reg[index]=A[j]
    Z.reg[index]=Z[[j]]
    ## Adjusted D and Y: D_ij* = n_jk * J * D_ij / n_k, Y_ij* similarly
    for (i in seq_along(index)){
      kk=B[[j]][i]
      weight_adj=njk[j,kk]*J/nk[kk]
      D.reg[index[i]]=D[[j]][i]*weight_adj
      Y.reg[index[i]]=Y[[j]][i]*weight_adj
    }
    ## Centered strata indicators: pi_k = I(B_ij=k) - n_k/N
    for (kk in 1:K){
      Pi.full[index, kk]=as.numeric(B[[j]]==kk) - nk[kk]/N.total
    }
  }
  ## Drop last stratum to avoid collinearity (sum of Pi columns = 0)
  Pi=Pi.full[,2:K,drop=FALSE]
  
  X= cbind(A.reg, 1-A.reg,  Z.reg*A.reg, Z.reg*(1-A.reg),
           Pi*A.reg, Pi*(1-A.reg))

  reg1s=lm(D.reg~0+X,weights=W)
  D.hat=as.vector(fitted(reg1s))

  M= cbind(A.reg, 1-A.reg,  D.hat*A.reg, D.hat*(1-A.reg),
           Pi*A.reg, Pi*(1-A.reg))
  reg2s=lm(Y.reg~0+M,weights=as.vector(W))
  res= Y.reg-cbind(A.reg, 1-A.reg,  D.reg*A.reg, D.reg*(1-A.reg),
                    Pi*A.reg, Pi*(1-A.reg))%*%reg2s$coefficients
  ###  variance
  p=4+2*(K-1)

  ## cluster robust variance
  MM=t(M)%*%diag(W)%*%M

  var.cluster.med=array(0,dim=c(p,p))
  for( j in 1:J){
    index= index.l[j]:index.r[j]
    Mj= M[index,]
    if (A[j]==1){
      Sj=cbind(W[index],   0, W[index]*D.hat[index], 0,
               W[index]*Pi[index,,drop=FALSE], matrix(0,nrow=length(index),ncol=K-1))

      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
    }else{
      Sj=cbind(0, W[index],   0, W[index]*D.hat[index],
               matrix(0,nrow=length(index),ncol=K-1), W[index]*Pi[index,,drop=FALSE])

      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])

    }
  }
  var.cluster=solve(MM)%*%var.cluster.med%*%solve(MM)


  ## cluster robust hc2 variance
  MM=t(M)%*%diag(W)%*%M

  var.cluster.med=array(0,dim=c(p,p))
  for( j in 1:J){
    index= index.l[j]:index.r[j]
    Mj= M[index,]
    if (A[j]==1){
      Sj=cbind(W[index],   0, W[index]*D.hat[index], 0,
               W[index]*Pi[index,,drop=FALSE], matrix(0,nrow=length(index),ncol=K-1))*sqrt(J1/(J1-1))

      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
    }else{
      Sj=cbind(0, W[index],   0, W[index]*D.hat[index],
               matrix(0,nrow=length(index),ncol=K-1), W[index]*Pi[index,,drop=FALSE])*sqrt((J-J1)/(J-J1-1))

      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])

    }
  }

  var.cluster.hc2=solve(MM)%*%var.cluster.med%*%solve(MM)
  ### individual robust hc2
  res.ind=rep(0,sum(n))
  var.ind.med=array(0,dim=c(p,p))
  for (j in 1:J){
    index= index.l[j]:index.r[j]
    adj1=   sum(res[index]*Z[[j]]/sum(Z[[j]]))
    adj0=   sum(res[index]*(1-Z[[j]])/sum(1-Z[[j]]))

    res.ind[index]=res[index] - ifelse(Z[[j]]==1,adj1,adj0)
  }

  for (j in 1:J){
    for (i in 1:n[j]){
      index=index.l[j]-1+i
      var.ind.med=var.ind.med+(M[index,])%*% t( M[index,])  *W[index]^2 * ifelse(Z.reg[index]==1, n1[j]/(n1[j]-1),(n0[j])/(n0[j]-1))*res.ind[index]^2
    }
  }
  var.ind=solve(MM)%*%var.ind.med%*%solve(MM)


  ### traditional hc2 variance
  var.hc2.med=array(0,dim=c(p,p))
  for (j in 1:J){
    for (i in 1:n[j]){
      index=index.l[j]-1+i
      if (A[j]==1){
        constant=ifelse(Z.reg[index]==1, J1*n1[j]/(J1*n1[j]-1),J1*n0[j]/(J1*n0[j]-1))
      }else{
        constant=ifelse(Z.reg[index]==1, J0*n1[j]/(J1*n1[j]-1),J0*n0[j]/(J1*n0[j]-1))
      }

      var.hc2.med=var.hc2.med+(M[index,])%*% t( M[index,])  *W[index]^2 * constant*res[index]^2
    }
  }
  var.hc2=solve(MM)%*%var.hc2.med%*%solve(MM)
  
  ## results
  est.CADE1=reg2s$coefficients[3]
  est.CADE0=reg2s$coefficients[4]
  var1.cluster=var.cluster[3,3]
  var0.cluster=var.cluster[4,4]
  var1.cluster.hc2=var.cluster.hc2[3,3]
  var0.cluster.hc2=var.cluster.hc2[4,4]
  var1.ind=var.ind[3,3]
  var0.ind=var.ind[4,4]
  var1.reg=(1-J1/J)*var.cluster.hc2[3,3]+(J1/J)*var.ind[3,3]
  var0.reg=(J1/J)*var.cluster.hc2[4,4]+(1-J1/J)*var.ind[4,4]
  var1.hc2=var.hc2[3,3]
  var0.hc2=var.hc2[4,4]

  return(list(CADE1=est.CADE1,CADE0=est.CADE0, var1.clu=var1.cluster,var0.clu=var0.cluster,var1.clu.hc2=var1.cluster.hc2,var0.clu.hc2=var0.cluster.hc2,   var1.ind=var1.ind,var0.ind=var0.ind,var1.reg=var1.reg,var0.reg=var0.reg,var1.hc2=var1.hc2,var0.hc2=var0.hc2,
    reg1s=reg1s, reg2s=reg2s
  ))
}
r1 <- r(data)
e = c()
for (i in 1:100){
  data <- assign_strata(data)
  r1 <- r(data)
  e = c(e, r1$CADE0 - cade_a0)
}

mean(e)
sd(e)
sqrt(mean(e^2))
