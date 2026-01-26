library('mvtnorm')
library('tidyverse')

set.seed(202001003)

#### total variance sigma 

### four mechanism with proportions 20 40 60 80


GenData_po = function(n,theta,r,rho,sigma=1){
   J <- length(n)
   m <- dim(theta)[2]
   Yij <- array(0,dim=c(sum(n),2,m))
   Yj <- array(0,dim=c(J,2,m))
   
   ##  between-cluster and within-cluster variances
   sigmaB <-  sigma*r
   sigmaW <-  sigma*(1-r)
   
   
   ### cluster level
   for (a in 1:4){
   Yj[,1,a] <- rnorm(J, theta[1,a],  sqrt(sigmaB) )
   Yj[,2,a] <- rnorm(J, theta[2,a]+rho*(Yj[,1,a]-theta[1,a]),  sqrt((1-rho^2)*sigmaB) )
   }
   
   ### individual level
   cum.n <- c(0,cumsum(n))
   for ( j in 1:J){
     for (a in 1:m){
       binormal <- t(rmvnorm(n[j],c(Yj[j,1,a],Yj[j,2,a]),   array(c(1,rho,rho,1),dim=c(2,2))*sigmaW ) )
       
       Yij[ (cum.n[j]+1):cum.n[j+1],1,a] <- binormal[1,]
       Yij[ (cum.n[j]+1):cum.n[j+1],2,a] <- binormal[2,]
     }
     
   }
   
   ### finite population adjustment
   
   effect.adj <- apply(Yij,2:3,mean)-theta
   
   for (a in 1:m){
     for (z in 1:2){
       Yij[,z,a] <- Yij[,z,a]-effect.adj[z,a]
     }
   }


  return(Yij)
   
}


GenData_obs = function(n,Yij,pa,qa){
  J<- length(n)
  id <- rep(1:J,n)
  
  
  A.cluster <-  sample(rep(1:m, qa*J ),J)
  
  Z <- rep(0,sum(n))
  A <- rep(0,sum(n))
  Y <- rep(0,sum(n))
  cum.n <- c(0,cumsum(n))
  for (j in 1:J){
    Z[(cum.n[j]+1):cum.n[j+1]] <- sample(rep(c(1,0),c(n[j]*pa[A.cluster[j]],n[j]*(1-pa[A.cluster[j]])+1e-7 ) ),n[j])
    A[(cum.n[j]+1):cum.n[j+1]] <- A.cluster[j]
    for (i in 1:n[j]){
      Y[cum.n[j]+i] <- Yij[cum.n[j]+i,Z[cum.n[j]+i]+1,A.cluster[j]]
    }
    
  }
  
  data <- data.frame(id=id,Z=Z,A=A,Y=Y)
  return(data)
}





Sim_equal <-  function(n.avg,J, theta, r, rho,sigma=1,pa,qa,effect="DE",repetition=1000){
  
  m  <- length(pa)
  Ja <- qa*J
  ## test statistics
  T <- rep(0,repetition)
  
  n <- rep(n.avg,J)
  po <- GenData_po(n,theta,r,rho,sigma=1) 
  
  for (rep in 1:repetition){
    if(rep%%200 == 0){     
      print(paste(rep, "/",repetition, sep = ""))
    }
    
    ### generate  observed data
    data <- GenData_obs (n,po,pa,qa)
    rej <- Test2SRE(data$id,data$Z,data$A,data$Y,qa,effect=effect)
    
    T[rep] <- rej
    
  }
  
  power <- mean(T)
  return(power)
  
}

rcategorical=function(n, p)
{
  rv=runif(n, 0, 1)
  out=rep(1,n)
  l=length(p)-1
  for(i in 1:l)
  {
    out=ifelse(rv>sum(p[1:i]),i+1,out)
  }
  return(out)
}

Sim_unequal <-  function(n.avg,J, theta, r, rho,sigma=1,pa,qa,effect="DE",repetition=1000){
  m  <- length(pa)
  Ja <- qa*J
  ## test statistics
  T <- rep(0,repetition)
  
  ### generate unequal cluste size
  n.index <- rcategorical(J,c(0.25,0.1,0.1,0.1,0.2,0.25))
  n <- n.avg* c(0.5,0.75,1,1.5,2,2.5)[n.index]
  po <- GenData_po(n,theta,r,rho,sigma=1) 
  
  for (rep in 1:repetition){
    if(rep%%200 == 0){     
      print(paste(rep, "/",repetition, sep = ""))
    }
    
    ### generate  observed data
    data <- GenData_obs (n,po,pa,qa)
    rej <- Test2SRE(data$id,data$Z,data$A,data$Y,qa,effect=effect)
    
    T[rep] <- rej
    
  }
  
  power <- mean(T)
  return(power)
  
}


roundm=function(x,m){
  max((floor(x/m)+1)*m,2*m)
}


source('analyze2SRE.R')
#### required sample size
m=4
qa=rep(1/m,m)
pa=c(0.2,0.4,0.6,0.8)
mu = 0.3

r.set=(0:20)/20
n.r=length(r.set)


## cluster size n=20     

J.DE.n20 <- rep(-1,n.r)
J.MDE.n20<- rep(-1,n.r)
J.SE.n20 <- rep(-1,n.r)
for ( i in 1:n.r){
  sam <- Calsamplesize(mu, 20, qa, pa, r=r.set[i],sigma=1,alpha=0.05, beta=0.2)
  J.DE.n20[i]=sam[1]
  J.MDE.n20[i]=sam[2]
  J.SE.n20[i]=sam[3]
}






## cluster size n=100     

J.DE.n100 <- rep(-1,n.r)
J.MDE.n100 <- rep(-1,n.r)
J.SE.n100 <- rep(-1,n.r)
for ( i in 1:n.r){
  sam <- Calsamplesize(mu, 100, qa, pa, r=r.set[i],sigma=1,alpha=0.05, beta=0.2)
  J.DE.n100 [i] <- sam[1]
  J.MDE.n100 [i] <- sam[2]
  J.SE.n100 [i] <- sam[3]
}



###### power simulation

theta.DE <- array(0,dim=c(2,4))
theta.MDE <- array(0,dim=c(2,4))
theta.SE <- array(0,dim=c(2,4))

theta.DE[1,] <-  runif(4,min=-0.3,max=0.3)
theta.DE[2,] <- theta.DE[1,]+0.3

theta.MDE[1,] <-  runif(4,min=-0.3,max=0.3)
theta.MDE[2,] <- theta.MDE[1,]+c(0.12,0.48,0.24,0.36)

theta.SE[1,1:3] <-  runif(3,min=-0.15,max=0.15)
theta.SE[2,1:3] <-  runif(3,min=-0.15,max=0.15)
theta.SE[1,4] <- 0.3+min(theta.SE[1,1:3])
theta.SE[2,4] <- 0.3+min(theta.SE[2,1:3])

#### n=20  equal cluster size  correct parameter (rho = 0)

power.DE.n20.equal.rho0=rep(-1,n.r)
power.MDE.n20.equal.rho0=rep(-1,n.r)
power.SE.n20.equal.rho0=rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n20.equal.rho0[i] <- Sim_equal(20,roundm(J.DE.n20[i],m), theta.DE, r=r.set[i], rho=0,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
 power.MDE.n20.equal.rho0[i]  <- Sim_equal(20,roundm(J.MDE.n20[i],m), theta.MDE,r=r.set[i], rho=0, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n20.equal.rho0[i] <-  Sim_equal(20,roundm(J.SE.n20[i],m), theta.SE, r=r.set[i], rho=0, sigma=1,pa,qa,effect="SE",repetition=1000)
}


#### n=20  unequal cluster size  correct parameter (rho = 0)
power.DE.n20.unequal.rho0=rep(-1,n.r)
power.MDE.n20.unequal.rho0=rep(-1,n.r)
power.SE.n20.unequal.rho0=rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n20.unequal.rho0[i] <- Sim_unequal(20,roundm(J.DE.n20[i],m), theta.DE, r=r.set[i], rho=0,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n20.unequal.rho0[i]  <- Sim_unequal(20,roundm(J.MDE.n20[i],m), theta.MDE,r=r.set[i], rho=0, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n20.unequal.rho0[i] <-  Sim_unequal(20,roundm(J.SE.n20[i],m), theta.SE, r=r.set[i], rho=0, sigma=1,pa,qa,effect="SE",repetition=1000)
}



#### n=100  equal cluster size  correct parameter (rho = 0)

power.DE.n100.equal.rho0 <- rep(-1,n.r)
power.MDE.n100.equal.rho0 <- rep(-1,n.r)
power.SE.n100.equal.rho0 <- rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n100.equal.rho0[i] <- Sim_equal(100,roundm(J.DE.n100[i],m), theta.DE, r=r.set[i], rho=0,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n100.equal.rho0[i]  <- Sim_equal(100,roundm(J.MDE.n100[i],m), theta.MDE,r=r.set[i], rho=0, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n100.equal.rho0[i] <-  Sim_equal(100,roundm(J.SE.n100[i],m), theta.SE, r=r.set[i], rho=0, sigma=1,pa,qa,effect="SE",repetition=1000)
}



#### n=100  unequal cluster size  correct parameter (rho = 0)
power.DE.n100.unequal.rho0=rep(-1,n.r)
power.MDE.n100.unequal.rho0=rep(-1,n.r)
power.SE.n100.unequal.rho0=rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n100.unequal.rho0[i] <- Sim_unequal(100,roundm(J.DE.n100[i],m), theta.DE, r=r.set[i], rho=0,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n100.unequal.rho0[i]  <- Sim_unequal(100,roundm(J.MDE.n100[i],m), theta.MDE,r=r.set[i], rho=0, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n100.unequal.rho0[i] <-  Sim_unequal(100,roundm(J.SE.n100[i],m), theta.SE, r=r.set[i], rho=0, sigma=1,pa,qa,effect="SE",repetition=1000)
}




############  mis specified parameter 
#### n=20  equal cluster size  rho=0.3

power.DE.n20.equal.rho03=rep(-1,n.r)
power.MDE.n20.equal.rho03=rep(-1,n.r)
power.SE.n20.equal.rho03=rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n20.equal.rho03[i] <- Sim_equal(20,roundm(J.DE.n20[i],m), theta.DE, r=r.set[i], rho=0.3,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n20.equal.rho03[i]  <- Sim_equal(20,roundm(J.MDE.n20[i],m), theta.MDE,r=r.set[i], rho=0.3, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n20.equal.rho03[i] <-  Sim_equal(20,roundm(J.SE.n20[i],m), theta.SE, r=r.set[i], rho=0.3, sigma=1,pa,qa,effect="SE",repetition=1000)
}


#### n=20  unequal cluster size  rho=0.3
power.DE.n20.unequal.rho03=rep(-1,n.r)
power.MDE.n20.unequal.rho03=rep(-1,n.r)
power.SE.n20.unequal.rho03=rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n20.unequal.rho03[i] <- Sim_unequal(20,roundm(J.DE.n20[i],m), theta.DE, r=r.set[i], rho=0.3,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n20.unequal.rho03[i]  <- Sim_unequal(20,roundm(J.MDE.n20[i],m), theta.MDE,r=r.set[i], rho=0.3, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n20.unequal.rho03[i] <-  Sim_unequal(20,roundm(J.SE.n20[i],m), theta.SE, r=r.set[i], rho=0.3, sigma=1,pa,qa,effect="SE",repetition=1000)
}


#### n=20  equal cluster size  rho=0.6

power.DE.n20.equal.rho06=rep(-1,n.r)
power.MDE.n20.equal.rho06=rep(-1,n.r)
power.SE.n20.equal.rho06=rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n20.equal.rho06[i] <- Sim_equal(20,roundm(J.DE.n20[i],m), theta.DE, r=r.set[i], rho=0.6,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n20.equal.rho06[i]  <- Sim_equal(20,roundm(J.MDE.n20[i],m), theta.MDE,r=r.set[i], rho=0.6, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n20.equal.rho06[i] <-  Sim_equal(20,roundm(J.SE.n20[i],m), theta.SE, r=r.set[i], rho=0.6, sigma=1,pa,qa,effect="SE",repetition=1000)
}


#### n=20  unequal cluster size  correct parameter rho=0.6
power.DE.n20.unequal.rho06=rep(-1,n.r)
power.MDE.n20.unequal.rho06=rep(-1,n.r)
power.SE.n20.unequal.rho06=rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n20.unequal.rho06[i] <- Sim_unequal(20,roundm(J.DE.n20[i],m), theta.DE, r=r.set[i], rho=0.6,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n20.unequal.rho06[i]  <- Sim_unequal(20,roundm(J.MDE.n20[i],m), theta.MDE,r=r.set[i], rho=0.6, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n20.unequal.rho06[i] <-  Sim_unequal(20,roundm(J.SE.n20[i],m), theta.SE, r=r.set[i], rho=0.6, sigma=1,pa,qa,effect="SE",repetition=1000)
}


#### n=100  equal cluster size  rho=0.3

power.DE.n100.equal.rho03 <- rep(-1,n.r)
power.MDE.n100.equal.rho03 <- rep(-1,n.r)
power.SE.n100.equal.rho03 <- rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n100.equal.rho03[i] <- Sim_equal(100,roundm(J.DE.n100[i],m), theta.DE, r=r.set[i], rho=0.3,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n100.equal.rho03[i]  <- Sim_equal(100,roundm(J.MDE.n100[i],m), theta.MDE,r=r.set[i], rho=0.3, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n100.equal.rho03[i] <-  Sim_equal(100,roundm(J.SE.n100[i],m), theta.SE, r=r.set[i], rho=0.3, sigma=1,pa,qa,effect="SE",repetition=1000)
}



#### n=100  unequal cluster size  rho=0.3
power.DE.n100.unequal.rho03 <- rep(-1,n.r)
power.MDE.n100.unequal.rho03 <- rep(-1,n.r)
power.SE.n100.unequal.rho03 <- rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n100.unequal.rho03[i] <- Sim_unequal(100,roundm(J.DE.n100[i],m), theta.DE, r=r.set[i], rho=0.3,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n100.unequal.rho03[i]  <- Sim_unequal(100,roundm(J.MDE.n100[i],m), theta.MDE,r=r.set[i], rho=0.3, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n100.unequal.rho03[i] <-  Sim_unequal(100,roundm(J.SE.n100[i],m), theta.SE, r=r.set[i], rho=0.3, sigma=1,pa,qa,effect="SE",repetition=1000)
}


#### n=100  equal cluster size  rho=0.6

power.DE.n100.equal.rho06 <- rep(-1,n.r)
power.MDE.n100.equal.rho06 <- rep(-1,n.r)
power.SE.n100.equal.rho06 <- rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n100.equal.rho06[i] <- Sim_equal(100,roundm(J.DE.n100[i],m), theta.DE, r=r.set[i], rho=0.6,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n100.equal.rho06[i]  <- Sim_equal(100,roundm(J.MDE.n100[i],m), theta.MDE,r=r.set[i], rho=0.6, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n100.equal.rho06[i] <-  Sim_equal(100,roundm(J.SE.n100[i],m), theta.SE, r=r.set[i], rho=0.6, sigma=1,pa,qa,effect="SE",repetition=1000)
}



#### n=100  unequal cluster size  rho=0.6
power.DE.n100.unequal.rho06 <- rep(-1,n.r)
power.MDE.n100.unequal.rho06 <- rep(-1,n.r)
power.SE.n100.unequal.rho06 <- rep(-1,n.r)


for ( i in 1:n.r){
  print(paste("DE ",i, "/",21, sep = ""))
  power.DE.n100.unequal.rho06[i] <- Sim_unequal(100,roundm(J.DE.n100[i],m), theta.DE, r=r.set[i], rho=0.6,sigma=1,pa,qa,effect="DE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("MDE ",i, "/",21, sep = ""))
  power.MDE.n100.unequal.rho06[i]  <- Sim_unequal(100,roundm(J.MDE.n100[i],m), theta.MDE,r=r.set[i], rho=0.6, sigma=1,pa,qa,effect="MDE",repetition=1000)
}

for ( i in 1:n.r){
  print(paste("SE ",i, "/",21, sep = ""))
  
  power.SE.n100.unequal.rho06[i] <-  Sim_unequal(100,roundm(J.SE.n100[i],m), theta.SE, r=r.set[i], rho=0.6, sigma=1,pa,qa,effect="SE",repetition=1000)
}







save.image('/Users/Zhichao/Dropbox/github/kosuke/mismatch/interference2/code/simulation.RData')
load('/Users/Zhichao/Dropbox/problems/interference/code2/simulation_samplesize.RData')  

