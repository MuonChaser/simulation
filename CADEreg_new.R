#' Regression-based method for the complier average direct effect
#'
#'
#' This function computes the point estimates of the complier average direct effect (CADE) and four
#'  different variance estimates: the HC2 variance, the cluster-robust variance, the cluster-robust HC2
#'  variance and the variance proposed in the reference. The estimators calculated using this function
#'  are cluster-weighted, i.e., the weights are equal for each cluster. To obtain the indivudal-weighted
#'  estimators, please multiply the recieved treatment and the outcome by \code{n_jJ/N}, where
#'  \code{n_j} is the number of individuals in cluster \code{j}, \code{J} is the number of clusters and
#'  \code{N} is the total number of individuals.
#'
#'
#' For the details of the method implemented by this function, see the
#' references.
#'
#' @param data  A data frame containing the relevant variables. The names for the variables should be: ``Z'' for the treatment assignment,  ``D''  for the actual received treatment, ``Y'' for the outcome, ``A'' for the treatment assignment mechanism and ``id'' for the cluster ID. The variable for the cluster id should be a factor.
#' @return A list of class \code{CADEreg} which contains the following items:
#' \item{CADE1}{ The point estimate of CADE(1).  } \item{CADE0}{ The point estimate of CADE(0).  }
#' \item{var1.clu}{ The cluster-robust variance of CADE(1).   } \item{var0.clu}{ The cluster-robust variance of CADE(0).  }
#'\item{var1.clu.hc2}{ The cluster-robust HC2 variance of CADE(1).   }
#'\item{var0.clu.hc2}{ The cluster-robust HC2 variance of CADE(0).    }
#'\item{var1.hc2}{ The  HC2 variance of CADE(1).    }
#'\item{var0.hc2}{ The  HC2 variance of CADE(0).    }
#'\item{var1.ind}{ The  individual-robust variance of CADE(1).    }
#'\item{var0.ind}{ The  individual-robust variance of CADE(0).    }
#'\item{var1.reg}{ The  proposed variance of CADE(1).    }
#'\item{var0.reg}{ The  proposed variance of CADE(0).    }
#' @author Kosuke Imai, Department of Government and Department of Statistics, Harvard University
#' \email{imai@@Harvard.Edu}, \url{https://imai.fas.harvard.edu};
#' Zhichao Jiang, Department of Politics, Princeton University
#' \email{zhichaoj@@princeton.edu}.
#' @references Kosuke Imai, Zhichao Jiang and Anup Malani (2021).
#' \dQuote{Causal Inference with Interference and Noncompliance in the Two-Stage Randomized Experiments}, \emph{Journal of the American Statistical Association}.
#' @keywords two-stage randomized experiments
#' @export CADEreg




CADEreg_new=function(data){
  ## transform the data into list
  if(!is.factor(data$id)){stop('The cluster_id should be a factor variable.')}
  cluster.id=unique(data$id)
  n.cluster=length(cluster.id)
  Z=vector("list", n.cluster)
  D=vector("list", n.cluster)
  Y=vector("list", n.cluster)
  B=vector("list", n.cluster)
  A=rep(0,n.cluster)

  ## split once instead of O(N*J) comparisons per cluster
  cl_data <- split(data, data$id)
  for (i in seq_along(cluster.id)){
    cld <- cl_data[[as.character(cluster.id[i])]]
    Z[[i]]=as.numeric(cld$Z)
    D[[i]]=as.numeric(cld$D)
    Y[[i]]=cld$Y
    B[[i]]=as.numeric(cld$B)
    if (length(unique(cld$A))!=1){
      stop( paste0('The assignment mechanism in cluster ',i,' should be the same.'))
    }
    A[i]=cld$A[1]
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
  # njk[j,k]: number of individuals in stratum k within cluster j
  njk=matrix(0, nrow=J, ncol=K)
  ## tabulate is faster than double loop with sum(B[[j]]==k)
  for (j in 1:J){
    njk[j,] <- tabulate(B[[j]], nbins=K)
  }
  # nk[k]: total number of individuals in stratum k
  nk <- colSums(njk)
  N.total=sum(n)

  ## cumsum instead of O(J^2) partial sums
  index.r <- cumsum(n)
  index.l <- c(1L, head(index.r, -1L) + 1L)

  W=rep(0, N.total)
  for (j in 1:J){
    index=index.l[j]:index.r[j]
    W[index]=ifelse(A[j]==1, 1/J1,1/(J-J1)) *   ifelse(Z[[j]]==1, 1/n1[j],1/(n0[j]))
  }
  ## Design matrix in the first stage
  A.reg= rep(0,N.total)
  Z.reg=rep(0,N.total)
  D.reg=rep(0,N.total)
  Y.reg=rep(0,N.total)
  Pi.full=matrix(0, nrow=N.total, ncol=K)
  for (j in 1:J){
    index=index.l[j]:index.r[j]
    A.reg[index]=A[j]
    Z.reg[index]=Z[[j]]
    ## Adjusted D and Y: vectorized weight adjustment
    Bj <- B[[j]]
    weight_adj_vec <- njk[j, Bj] * J / nk[Bj]
    D.reg[index] <- D[[j]] * weight_adj_vec
    Y.reg[index] <- Y[[j]] * weight_adj_vec
    ## Centered strata indicators: vectorized via outer product
    Pi.full[index, ] <- outer(Bj, seq_len(K), "==") -
                        matrix(nk / N.total, nrow = length(index), ncol = K, byrow = TRUE)
  }
  ## Drop last stratum to avoid collinearity (sum of Pi columns = 0)
  Pi=Pi.full[,2:K,drop=FALSE]

  X= cbind(A.reg, 1-A.reg,  Z.reg*A.reg, Z.reg*(1-A.reg),
           Pi*A.reg, Pi*(1-A.reg),
           Pi*Z.reg*A.reg, Pi*Z.reg*(1-A.reg))

  reg1s=lm(D.reg~0+X,weights=W)
  D.hat=as.vector(fitted(reg1s))

  M= cbind(A.reg, 1-A.reg,  D.hat*A.reg, D.hat*(1-A.reg),
           Pi*A.reg, Pi*(1-A.reg),
           Pi*Z.reg*A.reg, Pi*Z.reg*(1-A.reg))
  reg2s=lm(Y.reg~0+M,weights=as.vector(W))
  res= Y.reg-cbind(A.reg, 1-A.reg,  D.reg*A.reg, D.reg*(1-A.reg),
                    Pi*A.reg, Pi*(1-A.reg),
                    Pi*Z.reg*A.reg, Pi*Z.reg*(1-A.reg))%*%reg2s$coefficients
  ###  variance
  p=4+4*(K-1)

  ## Compute MM and MM_inv once; use crossprod to avoid creating N×N diag matrix
  MM  <- crossprod(M * as.vector(W), M)
  MM_inv <- solve(MM)

  ## cluster robust variance
  var.cluster.med=array(0,dim=c(p,p))
  for( j in 1:J){
    index= index.l[j]:index.r[j]
    Mj= M[index,]
    if (A[j]==1){
      Sj=cbind(W[index],   0, W[index]*D.hat[index], 0,
               W[index]*Pi[index,,drop=FALSE], matrix(0,nrow=length(index),ncol=K-1),
               W[index]*Pi[index,,drop=FALSE]*Z.reg[index], matrix(0,nrow=length(index),ncol=K-1))

      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
    }else{
      Sj=cbind(0, W[index],   0, W[index]*D.hat[index],
               matrix(0,nrow=length(index),ncol=K-1), W[index]*Pi[index,,drop=FALSE],
               matrix(0,nrow=length(index),ncol=K-1), W[index]*Pi[index,,drop=FALSE]*Z.reg[index])

      var.cluster.med=var.cluster.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])

    }
  }
  var.cluster=MM_inv%*%var.cluster.med%*%MM_inv


  ## cluster robust hc2 variance
  var.cluster.hc2.med=array(0,dim=c(p,p))
  for( j in 1:J){
    index= index.l[j]:index.r[j]
    Mj= M[index,]
    if (A[j]==1){
      Sj=cbind(W[index],   0, W[index]*D.hat[index], 0,
               W[index]*Pi[index,,drop=FALSE], matrix(0,nrow=length(index),ncol=K-1),
               W[index]*Pi[index,,drop=FALSE]*Z.reg[index], matrix(0,nrow=length(index),ncol=K-1))*sqrt(J1/(J1-1))

      var.cluster.hc2.med=var.cluster.hc2.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])
    }else{
      Sj=cbind(0, W[index],   0, W[index]*D.hat[index],
               matrix(0,nrow=length(index),ncol=K-1), W[index]*Pi[index,,drop=FALSE],
               matrix(0,nrow=length(index),ncol=K-1), W[index]*Pi[index,,drop=FALSE]*Z.reg[index])*sqrt((J-J1)/(J-J1-1))

      var.cluster.hc2.med=var.cluster.hc2.med+t(Sj)%*%res[index] %*%t(t(Sj)%*%res[index])

    }
  }
  var.cluster.hc2=MM_inv%*%var.cluster.hc2.med%*%MM_inv

  ### individual robust hc2
  res.ind=rep(0,N.total)
  var.ind.med=array(0,dim=c(p,p))
  for (j in 1:J){
    index= index.l[j]:index.r[j]
    adj1=   sum(res[index]*Z[[j]]/sum(Z[[j]]))
    adj0=   sum(res[index]*(1-Z[[j]])/sum(1-Z[[j]]))
    res.ind[index]=res[index] - ifelse(Z[[j]]==1,adj1,adj0)
  }

  ## Vectorize individual-level outer-product loops using crossprod
  cluster_idx <- rep(seq_len(J), n)
  n1_ind <- n1[cluster_idx]
  n0_ind <- n0[cluster_idx]

  hc2_scale_ind <- ifelse(Z.reg == 1, n1_ind / (n1_ind - 1), n0_ind / (n0_ind - 1))
  w_ind <- as.vector(W)^2 * hc2_scale_ind * as.vector(res.ind)^2
  var.ind.med <- crossprod(M * sqrt(w_ind))
  var.ind=MM_inv%*%var.ind.med%*%MM_inv

  ### traditional hc2 variance
  A_ind    <- A[cluster_idx]
  Ja_ind   <- ifelse(A_ind == 1, J1, J0)
  nz_ind   <- ifelse(Z.reg == 1, n1_ind, n0_ind)
  constant_vec <- Ja_ind * nz_ind / (J1 * nz_ind - 1)
  w_hc2 <- as.vector(W)^2 * constant_vec * as.vector(res)^2
  var.hc2.med <- crossprod(M * sqrt(w_hc2))
  var.hc2=MM_inv%*%var.hc2.med%*%MM_inv

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
  list(
    CADE1        = est.CADE1,
    CADE0        = est.CADE0,
    var1.clu     = var1.cluster,
    var0.clu     = var0.cluster,
    var1.clu.hc2 = var1.cluster.hc2,
    var0.clu.hc2 = var0.cluster.hc2,
    var1.ind     = var1.ind,
    var0.ind     = var0.ind,
    var1.reg     = var1.reg,       # proposed variance estimator
    var0.reg     = var0.reg,
    var1.hc2     = var1.hc2,
    var0.hc2     = var0.hc2,
    reg1s        = reg1s,
    reg2s        = reg2s
  )
}


