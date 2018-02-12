# Run SSA
# s             vector containing the initial state  
# reactionlist: List of strings containing reactions like "l1R1+l1R2->m1P1@0.1" or "1S+1E->1SE@0.1" 
#               first pattern is used by default, second pattern could be used if legend containing 
#               species names is provided in the legend and rename is set to THRUE  
#legend:        (optional, default empty) vector containing species names
#rename:        (optional, default FALSE) should pattern 2 be renamed to pattern one?
#filename:      (optional, default FALSE) name of the file 
ssa=function(s,reactionlist,T,legend=c(),rename=FALSE,filename=FALSE){
  # rename pattern 2 to 1
  if (length(legend)==length(s) & rename) { 
    reactionlist = rename_reactions(reactionlist,legend)  
  }
  smat=matrix(s,nrow=1,ncol=length(s))
  tvec=vector("numeric",1)
  t=0
  #convert reactionlist into matrix and c vector
  rvector=convreact(reactionlist,length(s))
  # reaction matrix with one row for each reaction 
  # positive numbers are products, negative reactants 
  R=matrix(as.numeric(rvector[1:(length(s)*length(reactionlist))]),nrow = length(reactionlist),ncol = length(s))
  # vector with stochastic rate constants for each reaction
  c=as.numeric(rvector[length(R)+1:length(reactionlist)])
  while (t<T){
    a=propensity(s,R,c)
    if (dim(a)[1]>0) {
      t=t+rexp(1,sum(a[,2]))
      if (dim(a)[1]==1) {
        s=s+R[a[,1],]  
      } else {
        s=s+R[sample(a[,1],1,prob=(a[,2]/sum(a[,2]))),]
      }
    }
    else {
      t=100  
    }
    smat=rbind(s,smat)
    tvec=c(tvec,t)
  }  
  matplot(tvec, smat[nrow(smat):1,], type = "s")
  if (length(legend)==length(s) & length(s)<5) { 
    legend("topright", legend=legend,col=seq_len(ncol(smat)),cex=1,fill=seq_len(ncol(smat)))
  }
  if (filename!=FALSE) {
    smat = cbind(tvec, smat[nrow(smat):1,]) 
    rownames(smat)=NULL
    colnames(smat)=NULL
    write.csv(smat, file = filename)
  }
}

# Renames reactions from pattern two to pattern one
# reactionlist: vector of strings containing reactions in pattern two
# legend:       vector containing species names used in pattern two
# returns:      vector of strings containing reactions in pattern one
rename_reactions=function(reactionlist, legend) {
  for(reaction in 1:length(reactionlist)) {
    for(name in 1:length(legend)) {
      posiname=gregexpr(sprintf("[[:digit:]]+%s[[:punct:]]",legend[name]), reactionlist[reaction])
      while (posiname[[1]][1]>0){
        posiarrow=regexpr("->",reactionlist[reaction])[1]
        start=substr(reactionlist[reaction],1,(posiname[[1]][1])-1)
        numberofs=substr(reactionlist[reaction],(posiname[[1]][1]),(posiname[[1]][1]+attr(posiname[[1]],"match.length")[1]-2-nchar(legend[name])))
        end=substr(reactionlist[reaction],(posiname[[1]][1]+attr(posiname[[1]],"match.length")[1]-1),nchar(reactionlist[reaction]))
        if(posiname[[1]][1]>posiarrow) {
          reactionlist[reaction]=sprintf("%sm%sP%s%s",start,numberofs,name,end)
        }
        else {
          reactionlist[reaction]=sprintf("%sl%sR%s%s",start,numberofs,name,end)
        }
        posiname=gregexpr(sprintf("[[:digit:]]+%s[[:punct:]]",legend[name]), reactionlist[reaction])
      }
    }  
  }
  reactionlist
}

# Converts vector of strings with reactions into reaction matrix and rate vector
# r:        vector of strings containing reactions in pattern one
# n:        number of species 
# returns:  vector containing reaction matreix and rate vector
convreact=function(r,n) {
  rm=matrix(nrow=0,ncol=n)
  k=vector("numeric",0)
  for (i in 1:length(r)) {
    if (dim(rm)[1]>1) {
      rm=rm[nrow(rm):1,]
    }
    rm=rbind(vector("numeric",n),rm)
    if (dim(rm)[1]>1) {
      rm=rm[nrow(rm):1,]
    }
    atsplit=strsplit(r[i],"@")
    k=c(k,atsplit[[1]][2])
    for (rp in 1:2) {
      arrowsplit=strsplit(atsplit[[1]][1],"->")[[1]][rp]
      if (!is.na(arrowsplit)) {
        plussplit=strsplit(arrowsplit,"\\+")
        for (spec in 1:length(plussplit[[1]])) {
          count = regmatches(plussplit[[1]][spec], gregexpr("[[:digit:]]+", plussplit[[1]][spec]))
          if (rp==1) {
            rm[i,as.numeric(count[[1]][2])]=as.numeric(count[[1]][1])*-1
          }
          else {
            rm[i,as.numeric(count[[1]][2])]=as.numeric(count[[1]][1])
          }
        }
      }
    }
    
  }
  c(rm,k)
}

# Calculates propensity functions for given state 
# s:        current state 
# R:        reaction matreix
# c:        rate vector
# returns:  matrix with propensity function results  
propensity=function(s,R,c){
  #IDs of posible reactions (where state + reaction is not -1)
  rid = which(apply(t(s+t(R)), 1, function(x) {min(x)>-1}))
  a=matrix(0,nrow=length(rid),ncol=2)
  #if reactins are possible
  if (length(rid)>0) {
    for (i in 1:length(rid)) {
      #IDs of reagents in the reaction
      sid=which(R[rid[i],]<0)
      a[i,] = c(rid[i],prod(choose(s[sid],R[rid[i],sid]*-1))*c[rid[i]])
    }
  }
  a
}



#################################
### Examples                  ###
#################################

# simple SSA test with pattern one from slides
ssa_test=function()  {
  s=c(0,10,8,0)
  reactionlist=c("l1R3+l1R2->m1P4@0.02",
                 "l1R1->@0.5",
                 "l1R4->m1R1@0.8")
  ssa(s,reactionlist,100)  
}

# Michaelis-Menten model with pattern one
mime_test=function()  {
  #converting ODE rate constants into stochastic
  #k1=1e6
  #k1r=1e-4
  #k2=0.1
  #vol=1e-15
  #NAc=6.02214129e23
  #c1=k1/(NAc*vol)
  #c1r=k1r/(NAc*vol)
  #c2=k2
  # other constans could be 0.001, 0.02, 0.8
  s=c(200,100,0,0)
  reactionlist=c("l1R1+l1R2->m1P3@0.001660539",
                 "l1R3->m1P1+m1P2@1.660539e-13",
                 "l1R3->m1R2+m1P4@0.1")
  ssa(s,reactionlist,100,c("S","E","SE","P"))  
}

# Michaelis-Menten model with pattern two (renamed)
mime_renamed_test=function() {
  s=c(2000,1000,0,0)
  reactionlist=c("1S+1E->1SE@0.001660539",
                 "1SE->1E+1S@1.660539e-13",
                 "1SE->1E+1P@0.1")
  ssa(s,reactionlist,100,c("S","E","SE","P"),rename = TRUE)  
}

# Riboflavin model with pattern two with 9 reactions
ribo_test1=function() {
  s=c(2,1,0,2,1,0,2,1,0,0)*1000
  legend=c("R5P","RibA","R5PRibA","D2B4P","RibH","RibHD2B4P","D8RL","RibE","D8RLRibE","Riboflavin")
  reactionlist=c("1R5P+1RibA->1R5PRibA@0.001660539",
                 "1R5PRibA->1R5P+1RibA@1.660539e-13",
                 "1R5PRibA->1RibA+1D2B4P@0.1",
                 
                 "1D2B4P+1RibH->1RibHD2B4P@0.001660539",
                 "1RibHD2B4P->1D2B4P+1RibH@1.660539e-13",
                 "1RibHD2B4P->1RibH+1D8RL@0.1",
                 
                 "1D8RL+1RibE->1D8RLRibE@0.001660539",
                 "1D8RLRibE->1D8RL+1RibE@1.660539e-13",
                 "1D8RLRibE->1RibE+1Riboflavin@0.1"
  )
  ssa(s,reactionlist,100,legend,rename = TRUE)  
}

# Riboflavin model with pattern two with 9 reactions, saves them to file
ribo_testfile=function() {
  s=c(2,1,0,2,1,0,2,1,0,0)*100
  legend=c("R5P","RibA","R5PRibA","D2B4P","RibH","RibHD2B4P","D8RL","RibE","D8RLRibE","Riboflavin")
  reactionlist=c("1R5P+1RibA->1R5PRibA@0.001660539",
                 "1R5PRibA->1R5P+1RibA@1.660539e-13",
                 "1R5PRibA->1RibA+1D2B4P@0.1",
                 
                 "1D2B4P+1RibH->1RibHD2B4P@0.001660539",
                 "1RibHD2B4P->1D2B4P+1RibH@1.660539e-13",
                 "1RibHD2B4P->1RibH+1D8RL@0.1",
                 
                 "1D8RL+1RibE->1D8RLRibE@0.001660539",
                 "1D8RLRibE->1D8RL+1RibE@1.660539e-13",
                 "1D8RLRibE->1RibE+1Riboflavin@0.1")
  ssa(s,reactionlist,100,legend,rename = TRUE, filename = "Riboflavin.csv")  
}

# Riboflavin model with pattern two with 15 reactions
ribo_test2=function() {
  s=c(2,1,0,2,1,0,2,1,0,0,1,0,0,0,0)*100
  legend=c("R5P","RibA","R5PRibA","D2B4P","RibH","RibHD2B4P","D8RL",
           "RibE","D8RLRibE","Riboflavin","RibC","FMN","FAD","RiboflavinRibC","FMNRibC")
  reactionlist=c("1R5P+1RibA->1R5PRibA@0.001660539",
                 "1R5PRibA->1R5P+1RibA@1.660539e-13",
                 "1R5PRibA->1RibA+1D2B4P@0.1",
                 
                 "1D2B4P+1RibH->1RibHD2B4P@0.001660539",
                 "1RibHD2B4P->1D2B4P+1RibH@1.660539e-13",
                 "1RibHD2B4P->1RibH+1D8RL@0.1",
                 
                 "1D8RL+1RibE->1D8RLRibE@0.001660539",
                 "1D8RLRibE->1D8RL+1RibE@1.660539e-13",
                 "1D8RLRibE->1RibE+1Riboflavin@0.1",
                 
                 "1Riboflavin+1RibC->1RiboflavinRibC@0.001660539",
                 "1RiboflavinRibC->1Riboflavin+1RibC@1.660539e-13",
                 "1RiboflavinRibC->1RibC+1FMN@0.1",
                 
                 "1FMN+1RibC->1FMNRibC@0.001660539",
                 "1FMNRibC->1FMN+1RibC@1.660539e-13",
                 "1FMNRibC->1RibC+1FAD@0.1"
                 )
  ssa(s,reactionlist,100,legend,rename = TRUE)  
}
