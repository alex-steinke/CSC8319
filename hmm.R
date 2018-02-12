# Generates random Lambda of size n.
# n         number of transitions 
# range:    Upper range for transitions probabilities (20 for real lambda 100 for first estimate) 
# returns:  matrix of size n containing transition probabilities  
gen_lam=function(n,range){
  lambda=matrix(0,ncol=n,nrow=n)
  for (i in 1:n) {
    for (j in 1:n) {
      if (i==j) {
        lambda[i,j] = sample(1:10000,1)   
      }
      else {
        lambda[i,j] = sample(1:range,1)
      }
    }  
    lambda[i,] = lambda[i,]/sum(lambda[i,])
  }
  lambda
}

# Generates random emission probabilities for Lambda. 
# alph      Number of symbols 
# ord       HMM order (only 0 or 1 works)
# trans     number of states in Lambda
# returns:  3-dimensinal array with the structure [alph, alph^ord, trans] containing emission probabilities
gen_p=function(alph, ord, trans){
  p = array(rep(0, alph*alph^ord*trans), dim=c(alph, alph^ord, trans))
  for (t in 1:trans) {
    for (a in 1:alph) {
      for (o in 1:alph^ord) {
        p[a,o,t] = sample(1:20,1)  
      }
      if (ord!=0) {
        p[a,,t] = p[a,,t]/sum(p[a,,t])  
      }
    }
    if (ord==0) {
      p[,,t] = p[,,t]/sum(p[,,t]) 
    }
  }
  p
}

# Calculates equilibrium probabilities. 
# P         Matrix of size n*n
# returns:  vector with equilibrium probabilities of P
equil=function(P){
  e=eigen(t(P))$vectors[,1]
  # to avoid crash in global bacause max() does not support complex numbers
  if (is.complex(e)) {
    e=Re(e)
  }
  e/sum(e)
}

# Generate random hidden sequence.
# n         length of the sequence
# lambda:   transition probabilities 
# returns:  hidden sequence
hssim=function(n,lambda)
{
  r=dim(lambda)[1]
  states=1:r
  s=vector("numeric",n)
  pi.lam=equil(lambda)
  s[1]=sample(states,1,FALSE,pi.lam)
  for (t in 2:n) {
    s[t]=sample(states,1,FALSE,prob=lambda[s[t-1],])
  }
  s
}

# Generate observed sequence.
# s         hidden sequence
# P:        emission probabilities
# returns:  observed sequence
obssim=function(s,P)
{
  n=length(s)
  q=dim(P)[1]
  r=dim(P)[3]
  states=1:q
  obs=vector("numeric",n)
  # 1st order HMM
  if (dim(P)[2]>1) {
    obs[1] = sample(states,1,FALSE,prob=equil(P[,,s[1]]))
    for (t in 2:n) {
      obs[t]=sample(states,1,FALSE,prob=P[,obs[t-1],s[t]])
    }
  }
  # 0 order HMM
  else {
    for (t in 1:n) {
      obs[t]=sample(states,1,FALSE,prob=P[,1,s[t]])
    }
  }
  obs
}

# Calculates the probability of observed symbol.
# P:        emission probabilities
# obs:      observed sequence
# obspos:   current position in the sequence
# returns:  vector of probabilities for observed symbol in all states
forwards_helper=function(P, obs, obspos){
  r=vector("numeric",dim(P)[3])
  # 1st order HMM
  if (dim(P)[2]>1) {
    # for 1st position use equilibrium probabilities
    if (obspos==1) {
      for (i in 1:dim(P)[3]) {
        r[i] = equil(P[,,i])[obs[1]]
      }
    }  
    else {
      r = P[obs[obspos],obs[obspos-1],]
    }
  }
  # 0 order HMM
  else {
    r=P[obs[obspos],1,]
  }
  r
}

# Estimate the hidden sequence using local decoding. 
# obs:      observed sequence
# lambda:   transition probabilities 
# P:        emission probabilities
# returns:  estimated hidden sequence
local=function(obs,lambda,P)
{
  r=dim(lambda)[1]
  q=dim(P)[1]
  n=length(obs)
  s=vector("numeric",n)
  f=matrix(0,nrow=r,ncol=n)
  b=matrix(0,nrow=r,ncol=n)
  # forwards
  f0=equil(lambda)
  f[,1]=forwards_helper(P,obs,1)*(f0%*%lambda)
  for (i in 2:n) {
    f[,i]=forwards_helper(P,obs,i)*(f[,i-1]%*%lambda)
    f[,i]=f[,i]/sum(f[,i])
  }
  # backwards
  b[,n]=rep(1,r)
  s[n]=which.max(f[,n]*b[,n])
  for (i in (n-1):1) {
    if (dim(P)[2]>1) {
      b[,i]=(P[obs[i+1],obs[i],]*b[,i+1])%*%lambda
    }
    else {
      b[,i]=(P[obs[i+1],1,]*b[,i+1])%*%lambda
    }
    b[,i]=b[,i]/sum(b[,i])
    s[i]=which.max(f[,i]*b[,i])
  }
  s
}

# Estimate the hidden sequence using global decoding. 
# obs:      observed sequence
# lambda:   transition probabilities 
# P:        emission probabilities
# returns:  estimated hidden sequence
global=function(obs,lambda,P)
{
  r=dim(lambda)[1]
  q=dim(P)[1]
  n=length(obs)
  s=vector("numeric",n)
  f=matrix(0,nrow=r,ncol=n)
  b=matrix(0,nrow=r,ncol=n)
  # forwards
  f0=equil(lambda)
  f[,1]=forwards_helper(P,obs,1)*(f0%*%lambda)
  for (i in 2:n) {
    for (k in 1:r) {
      f[k,i]=forwards_helper(P,obs,i)[k]*max(f[,i-1]*lambda[,k])
    }
    f[,i]=f[,i]/sum(f[,i])
  }	
  # backwards
  s[n]=which.max(f[,n])
  for (i in (n-1):1) {
    s[i]=which.max(lambda[,s[i+1]]*f[,i])
  }
  s
}


# Estimate parameters using changed version of local decoding
# obs:      observed sequence
# lambda:   transition probabilities 
# P:        emission probabilities
# returns:  vector containing new lambda and P
paramest=function(obs,lambda,P){
  r=dim(lambda)[1]
  q=dim(P)[1]
  q2=dim(P)[2]
  n=length(obs)
  s=vector("numeric",n)
  f=matrix(0,nrow=r,ncol=n)
  m=matrix(0,nrow=r,ncol=n)
  b=array(rep(0, r,r,n), dim=c(r,r,n))
  lambda.hat = lambda
  lambda.hat[,]=0
  p.hat = P
  p.hat[,,]=0
  # temporary sums for division
  lam.div=rep(0,r)
  p.div=matrix(0,nrow=q,ncol=r)
  f0=equil(lambda)
  f[,1]=forwards_helper(P,obs,1)*(f0%*%lambda)
  for (i in 2:n) {
    f[,i]=forwards_helper(P,obs,i)*(f[,i-1]%*%lambda)
    f[,i]=f[,i]/sum(f[,i])
  }
  # calculate b[j,k,t] and m[j,t] using backwards part
  m[,n]=f[,n]
  for (t in (n):1) {
    for (j in 1:r) {
      for (k in 1:r) {
        b[j,k,t]=lambda[j,k]*f[j,t]/sum(lambda[,k]*f[,t])
      } 
      if (t<n) {
        m[j,t]=b[j,,t]%*%m[,t+1]    
      }
    }
    
  }
  #calculate lambda.hat and p.hat using b[j,k,t] and m[j,t]
  for (j in 1:r) {
    for (t in 2:n) {
      for (k in 1:r) {
        lambda.hat[j,k]=lambda.hat[j,k]+b[j,k,t-1]*m[k,t]
        p.hat[obs[t],obs[t-1],k]=p.hat[obs[t],obs[t-1],k]+m[k,t]  
        p.div[obs[t],k]=p.div[obs[t],k]+m[k,t]  
      } 
      lam.div[j]=lam.div[j]+m[j,t-1]
    }
  }
  lambda.hat=lambda.hat/lam.div
  for (j in 1:q) {
    for (k in 1:r) {
      p.hat[j,,k]=p.hat[j,,k]/p.div[j,k]  
    }  
  }
  c(lambda.hat,p.hat)
}

# Calculate the accuracy of estimations.
# h:        hidden sequence
# l:        estimated hidden sequence using local decoding
# g:        estimated hidden sequence using global decoding
my_count=function(h,l,g=0)
{
  counter_g=0
  counter_l=0
  for (i in 1:length(h)) {
    if(length(g)>1) {
      if (g[i]!=h[i]) {
        counter_g=counter_g+1  
      }  
    }
    if (l[i]!=h[i]) {
      counter_l=counter_l+1  
    }
  }
  print(paste("Errors local decoding: ",counter_l/10000*100," %"))
  if(length(g)>1) {
    print(paste("Errors global decoding: ",counter_g/10000*100," %"))
  }
}

# Replaces given P values with more specific for up to 5 states  
# Just for testing purposes   
fixp=function(P) {
  p.old=P
  P = array(rep(0, 4,4,5), dim=c(4,4,5))
  P[1,,1] = c(0.6,0.1,0.1,0.2)
  P[2,,1] = c(0.1,0.7,0.1,0.1)
  P[3,,1] = c(0.1,0.1,0.5,0.3)
  P[4,,1] = c(0.2,0.3,0.1,0.4)
  
  P[1,,2] = c(0.1,0.1,0.1,0.7)
  P[2,,2] = c(0.2,0.1,0.6,0.2)
  P[3,,2] = c(0.1,0.5,0.1,0.3)
  P[4,,2] = c(0.7,0.1,0.1,0.1)
  
  P[1,,3] = c(0.1,0.4,0.4,0.1)
  P[2,,3] = c(0.3,0.4,0.1,0.2)
  P[3,,3] = c(0.6,0.1,0.2,0.1)
  P[4,,3] = c(0.1,0.1,0.7,0.1)
  
  P[1,,4] = c(0.3,0.3,0.3,0.1)
  P[2,,4] = c(0.1,0.3,0.3,0.3)
  P[3,,4] = c(0.3,0.3,0.1,0.3)
  P[4,,4] = c(0.3,0.3,0.3,0.1)
  
  P[1,,5] = c(0.1,0.5,0.3,0.1)
  P[2,,5] = c(0.1,0.3,0.5,0.1)
  P[3,,5] = c(0.1,0.3,0.5,0.1)
  P[4,,5] = c(0.1,0.5,0.3,0.1)
  for (i in 1:dim(p.old)[3]) {
    p.old[,,i]=P[,,i]
  }
  p.old
}

# Rename transitions in estimated sequence according to probabilities in hidden sequence 
# sh:       hidden sequence
# sh.est:   estimated hidden sequence
# returns:  estimated hidden sequence with renamed transitions
rename_est=function(sh, sh.est) {
  sh_count = as.data.frame(table(sh))
  sh_count[,1]=seq(1:length(sh_count[,2]))
  sh_count=sh_count[order(sh_count[,2],decreasing=TRUE),]
  shest_count = as.data.frame(table(sh.est))
  shest_count[,1]=seq(1:length(shest_count[,2]))
  shest_count=shest_count[order(shest_count[,2],decreasing=TRUE),]
  for (i in 1:length(shest_count[,2])) {
    sh.est[sh.est==shest_count[i,1]]=sh_count[i,1]*-1
  }
  sh.est=sh.est*-1
}

# Measures convergence  using RMS
# lam:      old lambda
# p:        old emission probabilities 
# lam.est:  estimated lambda
# p.est:    estimated emission probabilities 
# returns:  convergance
rms=function(lam,p,lam.est,p.est){
  old = c(p,lam)
  new = c(p.est,lam.est)
  sqrt(sum((new-old)^2)/length(old))
}

# Measures convergence  using MAD
# lam:      old lambda
# p:        old emission probabilities 
# lam.est:  estimated lambda
# p.est:    estimated emission probabilities 
# returns:  convergance
mad=function(lam,p,lam.est,p.est){
  old = c(p,lam)
  new = c(p.est,lam.est)
  sum(abs(new-old))/length(old)
}

showgraph=function(s,e,g=0){
  op=par(mfrow=c(3,1))
  plot(ts(s),main="Truth")
  plot(ts(e),main="Local decoding")
  if(length(g)>1) {
    plot(ts(g),main="Global decoding")
  }
  par(op)
}
# Estimates the hidden sequence using EM algorithm
# obs:      observed sequence
# measure:  convergence measure
em=function(obs, measure) {
  lam.est=rbind(gen_lam(trans_num,100))
  p.est=gen_p(alphabet_length,ord,trans_num)
  showgraph(sh,locest)
  conv = 1
  counter=0
  while (conv>10^-6 & counter<100 ) {
    counter=counter+1
    lam.old=lam.est
    p.old=p.est
    theta=paramest(obs,lam.old,p.old)
    lam.est=matrix(theta[1:length(lam.old)],nrow=dim(lam.old)[1],ncol=dim(lam.old)[1])
    p.est=array(theta[length(lam.old)+1:length(theta)], dim=dim(p.old))
    if (measure=="MAD") {
      conv=mad(lam.old,p.old,lam.est,p.est)  
    }
    else {
      conv=rms(lam.old,p.old,lam.est,p.est)    
    }
    print(paste("Iteration: ",counter,"     Convergence measure (",measure,"): ", conv))
    locest=local(obs,lam.est,p.est)
    showgraph(sh,locest)
  }
  locest=local(obs,lam.est,p.est)
  showgraph(sh,locest)
  locest=rename_est(sh,locest)
  showgraph(sh,locest)
  my_count(sh,locest)
  
}

# initial parameters 
seq_length = 10000
alphabet_length = 4
trans_num = 3
# HMM order, 0 or 1
ord = 1

# generate observed sequence
lambda=rbind(gen_lam(trans_num,20))
sh=hssim(seq_length,lambda)
P=gen_p(alphabet_length,ord,trans_num)
#Overwrite random P with premade one for testing purposes    
P=fixp(P)
obs=obssim(sh,P)

print("Decoding with known P and Lambda")
locest=local(obs,lambda,P)
globest=global(obs,lambda,P)
my_count(sh,locest,globest)
showgraph(sh,locest,globest)

print("Local decoding with EM")
em(obs, "RMS")