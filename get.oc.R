
###########################################################################################
#  This file contains get.oc()  function to generate operating characteristics for        #
#  several interval designs including the CCD, mTPI, BOIN, Keyboard and UMPBI designs.    #
#                                                                                         #
#  The function was adapted based on the codes written by Suyu Liu and Ying Yuan.         #
#                                                                                         #
#  If have any questions, please contact Ruitao Lin (ruitaolin@gmail.com)                 #
#                                                                                         #
###########################################################################################
#
#                                        #############
#                                        # Arguments:#
#                                        #############
###########################################################################################
#     target:  the target toxicity rate
#     p.true:  the true toxicity rate for each dose
#     ncohort:   the total number of cohorts
#     cohortsize:   cohort size
#     startdose:  the starting dose of the trial
#     design:  1 is the CCD design, 2 is the mTPI design, 3 is the BOIN design, 
#              4 is the Keyboard design (or mTPI2 design), 5 is the UMPBI design
#              (Default values for design parameters of these designs are utilized)
#			cutoff.eli:  cutoff to eliminate the overly toxic dose for safety monitoring
#     ntrial: number of simulated trials
############################################################################################


get.oc <- function(target, p.true, ncohort, cohortsize, startdose=1, design=5, cutoff.eli=0.95, ntrial=1000)
{	 	
  ### simple error checking
  if(sum(seq(1,5)==design)==0) 
  {cat("Please specify the correct design! (1--CCD, 2--mTPI, 3--BOIN, 4--Keyboard, 5--UMPBI) \n");  return(1);}

  ####################################################################################################
  ## to make get.oc as self-contained function, we copied functions get.boundary() and select.mtd() here.
  get.boundary <- function(target, ncohort, cohortsize=3,design=1,cutoff.eli=0.95)
  {
    density1 <- function(p, n, m1, m2) {pbinom(m1, n, p)+1-pbinom(m2-1, n, p);}
    density2 <- function(p, n, m1) {1-pbinom(m1, n, p)};
    density3 <- function(p, n, m2) {pbinom(m2-1, n, p);}
    df <- function(p,gamma,n,target){
      l<-(n/(1-p))/(log(p/(1-p))-log(target/(1-target)))
      l<-l-(log(gamma)-n*log((1-p)/(1-target)))/p/(1-p)/(log(p/(1-p))-log(target/(1-target)))^2
    }
    
    ### simple error checking
    if(sum(seq(1,5)==design)==0) 
    {cat("Please specify the correct design! (1--CCD, 2--mTPI, 3--BOIN, 4--Keyboard, 5--UMPBI) \n");  return(1);}
    
    ### default parameter setting
    if (design==1){
      if (target<0.3) {lambda1=target-0.09; lambda2=target+0.09}
      else if (target<0.4) {lambda1=target-0.1; lambda2=target+0.1}
      else if (target<0.45){lambda1=target-0.12; lambda2=target+0.12}
      else {lambda1=target-0.13; lambda2=target+0.13}
    }
    
    if (design==2 || design==4){
      p.saf = target-0.05; 
      p.tox = target+0.05;}
    if (design==3){
      p.saf = target*0.6; 
      p.tox = target*1.4;
      lambda1  = log((1-p.saf)/(1-target))/log(target*(1-p.saf)/(p.saf*(1-target)));
      lambda2  = log((1-target)/(1-p.tox))/log(p.tox*(1-target)/(target*(1-p.tox)));}
    if (design==5) {c<-log(1.1)/3}
    
    ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
    npts = ncohort*cohortsize;
    ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
    for(n in (1:ncohort)*cohortsize)
    {
      ntrt = c(ntrt, n);
      if (design==1||design==3){
        cutoff1 = floor(n*lambda1)
        cutoff2 = floor(n*lambda2)+1
      } 
      
      if (design==5){
        gamma = exp(c*sqrt(n))
        p.tox = uniroot(df,c(target+0.0001,target*2),gamma=gamma,n=n,target=target)$root
        p.saf = uniroot(df,c(0.0001,target-0.0001),gamma=gamma,n=n,target=target)$root
        lambda1  = (log((1-p.saf)/(1-target))-log(gamma)/n)/log(target*(1-p.saf)/(p.saf*(1-target)));
        lambda2  = (log((1-target)/(1-p.tox))+log(gamma)/n)/log(p.tox*(1-target)/(target*(1-p.tox)));
        cutoff1 = floor(n*lambda1)
        cutoff2 = floor(n*lambda2)+1
      }
      
      if (design==2 || design==4){
        error.min=3;
        for (m1 in 0:floor(target*n)){
          for (m2 in ceiling(target*n):n){
            if (design==2){
              error1 = integrate(density1, lower=p.saf, upper=p.tox, n, m1, m2)$value/(p.tox-p.saf);
              error2 = integrate(density2, lower=0, upper=p.saf, n, m1)$value/p.saf;
              error3 = integrate(density3, lower=p.tox, upper=1, n, m2)$value/(1-p.tox);
            }
            if (design==4){
              epsilon=p.tox-p.saf
              error1 = integrate(density1, lower=p.saf, upper=p.tox, n, m1, m2)$value/epsilon;
              error2 = integrate(density2, lower=p.saf-epsilon, upper=p.saf, n, m1)$value/epsilon;
              error3 = integrate(density3, lower=p.tox, upper=p.tox+epsilon, n, m2)$value/epsilon;
            }
            error=error1+error2+error3;
            if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
          }
        }
      }
      b.e = c(b.e, cutoff1);
      b.d = c(b.d, cutoff2);
      
      elimineed=0; # indicating whether elimination is needed
      if(n<3) { elim = c(elim, NA); }  # require treating at least 3 patients before eliminating a dose
      else
      {
        for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
        {
          if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
        }
        if(elimineed==1) { elim = c(elim, ntox); }
        else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
      }
    }
    for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
    boundaries = rbind(ntrt, elim, b.d, b.e);
    rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
                             "Deescalate if # of DLT >=",  "Escalate if # of DLT <=");
    colnames(boundaries) = rep("", ncohort);
    
    return(boundaries);
  }
  
  
  select.mtd <- function(target, y, n, cutoff.eli=0.95)
  {
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    pava <- function (x, wt = rep(1, length(x))) 
    {
      n <- length(x)
      if (n <= 1) 
        return(x)
      if (any(is.na(x)) || any(is.na(wt))) {
        stop("Missing values in 'x' or 'wt' not allowed")
      }
      lvlsets <- (1:n)
      repeat {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) 
          break
        i <- min((1:(n - 1))[viol])
        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i + 1]
        ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
        x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
        lvlsets[ilvl] <- lvl1
      }
      x
    }
    ## determine whether the dose has been eliminated during the trial
    ndose=length(n);
    elimi=rep(0, ndose);
    for(i in 1:ndose)
    {
      if(n[i]>2) {if(1-pbeta(target, y[i]+1, n[i]-y[i]+1)>cutoff.eli) {elimi[i:ndose]=1; break;}}
    }
    
    if(elimi[1]==1) { selectdose=99; } ## no dose should be selected if the first dose is already very toxic
    else
    {
      nadmis = min(max(which(elimi==0)), max(which(n!=0))); ## the highest admissble (or un-eliminated) dose level
      ## poster mean and variance of toxicity probabilities using beta(0.005, 0.005) as the prior 
      phat = (y[1:nadmis]+0.005)/(n[1:nadmis]+0.01); 
      phat.var = (y[1:nadmis]+0.005)*(n[1:nadmis]-y[1:nadmis]+0.005)/((n[1:nadmis]+0.01)^2*(n[1:nadmis]+0.01+1))
      
      ## perform the isotonic transformation using PAVA
      phat = pava(phat, wt=1/phat.var) 
      phat = phat + (1:nadmis)*1E-10 ## break ties by adding an increasingly small number 
      selectdose = sort(abs(phat-target), index.return=T)$ix[1]  ## select dose closest to the target as the MTD
    }
    return(selectdose);  
  }
  
  ########################### end of subroutines  ########################################
  
  set.seed(6);
  ndose=length(p.true);	
  npts = ncohort*cohortsize;
  Y=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
  dselect = rep(0, ntrial); # store the selected dose level
  
  ## obtain dose escalation and deescalation boundaries
  temp=get.boundary(target, ncohort, cohortsize,design, cutoff.eli); 	
  b.e=temp[4,];   # escalation boundary
  b.d=temp[3,];   # deescalation boundary
  b.elim=temp[2,];  # elimination boundary
  
  ################## simulate trials ###################
  for(trial in 1:ntrial)
  {
    y<-rep(0, ndose);    ## the number of DLT at each dose level
    n<-rep(0, ndose);    ## the number of patients treated at each dose level
    earlystop=0;         ## indiate whether the trial terminates early
    d=startdose;         ## starting dose level
    elimi = rep(0, ndose);  ## indicate whether doses are eliminated
    
    for(i in 1:ncohort)  
    {  			
      ### generate toxicity outcome
      y[d] = y[d] + sum(runif(cohortsize)<p.true[d]);
      n[d] = n[d] + cohortsize;
      nc = n[d]/cohortsize;
      
      ## determine if the current dose should be eliminated
      if(!is.na(b.elim[nc]))
      {
        if(y[d]>=b.elim[nc]) 
        {      
          elimi[d:ndose]=1;
          if(d==1) {earlystop=1; break;} 
        }
      }
      
      ## dose escalation/de-escalation
      if(y[d]<=b.e[nc] && d!=ndose) { if(elimi[d+1]==0) d=d+1; }
      else if(y[d]>=b.d[nc] && d!=1) { d=d-1; }
      else { d=d; }
    }
    Y[trial,]=y;
    N[trial,]=n;
    if(earlystop==1) { dselect[trial]=99; }
    else  dselect[trial]=select.mtd(target, y, n, cutoff.eli);		
  }	
  
  # output results 
  selpercent=rep(0, ndose);
  for(i in 1:ndose) { selpercent[i]=sum(dselect==i)/ntrial*100; }
  cat("selection percentage at each dose level (%):\n");
  cat(formatC(selpercent, digits=1, format="f"), sep="  ", "\n");
  cat("number of patients treated at each dose level:\n");
  cat(formatC(apply(N,2,mean), digits=1, format="f"), sep ="  ", "\n");
  cat("number of toxicity observed at each dose level:\n");
  cat(formatC(apply(Y,2,mean), digits=1, format="f"), sep ="  ", "\n");
  cat("average number of toxicties:", formatC(sum(Y)/ntrial, digits=1, format="f"), "\n");
  cat("average number of patients:", formatC(sum(N)/ntrial, digits=1, format="f"), "\n");
}


