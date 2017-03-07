Var<-function(x) mean(x^2)-mean(x)^2 #scalar version
vVar<-function(x) apply(x^2,1,mean)-apply(x,1,mean)^2 #vector version

Cov<-function(x,y) mean(x*y)-mean(x)*mean(y)
vCov<-function(x,y) apply(t(x)*y,2,mean)-apply(x,1,mean)*mean(y) #x is a matrix and y is a vector

sqrtsum<-function(y) sqrt(sum(y^2))


# x1 is log2 frequencies for the 1st replicate of all timepts
# x2 is log2 frequencies for the 2nd replicate of all timepts
# ab1 is abundance thresholds for all 1st replicates
# ab2 is abundance thresholds for all 2nd replicates
fit_ac_fc<-function(x1,ab1,x2,ab2) { #badx is TRUE when x-value is bad
  
  er_ac<-1
  l<-0
  nx<-nrow(x1)
  
  # Prediction: useless1 and useless2 here should have the same values as bad1 and bad2 in the 
  # Number of Constructs Below Abundance Threshold section ... boxed code seems to calc same result
  # as code in that section, but by slightly diff (but equivalent) means
  # -------------------
  # t(x1) has timept as rows and constructs as columns.  
  # I can't for the life of me figure out why it is necessary to transpose x1, compare it to ab1, 
  # and then transpose it back ... since the ">ab1" evaluation will be at the individual cell level,
  # and the individual cells hold log2 frequency values no matter which way you set up the rows and cols,
  # and since no margin calculations are being done in the good1 calculation, .... WHY?
  good1<-t(t(x1)>ab1) # constructs as rows, timept for 1st replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2 freq for row/col combination is above relevant abundance threshold
  good2<-t(t(x2)>ab2) # constructs as rows, timept for 2nd replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2 freq for row/col combination is above relevant abundance threshold
  useless1<-apply(good1,1,sum)<2 # 1 = sum over rows--i.e., constructs.  Any construct that isn't above abundance threshold in 1st replicate in at least 2 timepoints has 1/TRUE in "useless" matrix
  useless2<-apply(good2,1,sum)<2
  
  # ab
  #print(sum(useless1))
  #print(sum(useless2))
  # PREDICTION CORRECT!
  # --------------------
  
  # I do not understand how the statements below actually *remove* anything from the good1 and good2 dfs ...
  # It seems like, for each constructs that isn't above abundance thresholds in at least two timepoints for this replicate, it sets that construct's "good" values to FALSE for *all* timepoints in this replicate
  good1[useless1,]<-FALSE #remove singletons
  good2[useless2,]<-FALSE #remove singletons
  
  # Note: here apply(goodX,1,sum) is NOT uselessX because goodX was changed above
  #allbad is true for all the constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments
  allbad<-apply(good1,1,sum)<2 & apply(good2,1,sum)<2 #in this case I have nothing to use in either experiment
  
  # nt = number of timepoints
  lambda1<-rep(0,nt)
  lambda2<-rep(0,nt)
  ac1<-x1[,1] #just a guess # log2 frequencies for all constructs for first timepoint in this replicate
  ac2<-x2[,1] #just a guess
  fc<-rep(0,nx) # nx = number of constructs
  
  # for 1 to number of constructs
  for (i in 1:nx) {
    # if this construct doesn't have at least two timepoints above abundance threshold in at least one replicate, ignore it and move on
    if (allbad[i]) next #from now on there is at least one good experiment
    
    # apparently "f" stands for fitness--which is covariance (see below)--and # stands for replicate
    # v stands for variance and # " "
    f1<-0
    f2<-0
    v1<-0
    v2<-0
    
    g1<-good1[i,] # g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the first replicate 
    if (sum(g1)>1) { #it's a good experiment # if there are at least two good timepoints for this construct in replicate 1
      mx1<-mean(x1[i,g1]) # get the mean of the log2 frequencies for all the good timepoints for this construct in replicate 1
      # time is GLOBAL vector of timepoints (numbers, usually in days, in ascending order)
      mt1<-mean(time[g1]) # get mean of timepoints (e.g. mean number of days) for all the good timepoints for this construct in replicate 1
      v1<-Var(time[g1]) # get variance of timepoints " "
      f1<-Cov(x1[i,g1],time[g1]) # f1 = covariance of log2 frequencies for all the good timepoints for this construct in replicate 1 with the timepoints for those good timepoints
    }
    
    # do the exact same thing as above, but for replicate 2
    g2<-good2[i,]
    if (sum(g2)>1) { #it's a good experiment
      mx2<-mean(x2[i,g2])
      mt2<-mean(time[g2])
      v2<-Var(time[g2])
      f2<-Cov(x2[i,g2],time[g2])
    }
    
    # fc[i] is the combined fitness (across replicates) for construct i
    fc[i]<-(f1+f2)/(v1+v2) #the combined fitness from replicate 1+2
    #fc remains defined up to an additive constant
    
    if (sum(g1)>1) { # if there are at least two good timepoints for this construct in replicate 1
      # ac1 = mean of log2 freqs for this construct for good timepoints for rep 1 - (mean of timepts for good timepoints for this construct for rep 1)*combined fitness across replicates for this construct
      ac1[i]<-mx1-fc[i]*mt1
    }
    
    # same as above but for replicate 2
    if (sum(g2)>1) {
      ac2[i]<-mx2-fc[i]*mt2
    }
  }
  
  # ac is the initial condition (in log2 frequency) for construct c 
  # this is normalizing ac1 ... I think this is what is going on:
  # Roman's methods say "By definition, log2 relative frequencies satisfy the constraint
  # sum over c of (2^xc) = 1 at all times"
  # ac1 is the set of log2 frequencies for all constructs in replicate 1 at "initial conditions", so 
  # ac1 must satisfy the above constraint.  IFF the above constraint is satisfied, -log2(1) = 0.
  # If the above constraint isn't satisfied, the value of -log2(sum over c of (2^ac)) is subtracted from
  # ac to ensure the constraint is satisfied.
  alpha<- -log2(sum(2^ac1)) 
  ac1<-ac1+alpha #enforce normalization at time=0, sum(2^ac)=1
  
  # same as above but for replicate 2
  alpha<- -log2(sum(2^ac2))
  ac2<-ac2+alpha #enforce normalization at time=0, sum(2^ac)=1
  
  # ok, the *expected* log2 frequency for construct x at time t is:
  # xc(t) = ac + fc*t - log2(sum over c of 2^(ac + fc*t))
  # apparently the below code calls the term starting with "-log2" "lambda".
  # Note that lambda is being calculated separately for each combination of timepoint+replicate
  
  # for 1 to number of timepoints
  for (i in 1:nt) {
    lambda1[i]<- -log2(sum(2^(ac1+fc*time[i])))
    lambda2[i]<- -log2(sum(2^(ac2+fc*time[i])))
  } #these are initial estimates of lambda(t)
  
  # x1 is log2 frequencies for the 1st replicate of all timepts
  xfit1<-x1 #for size # I think this means that xfitX is being set to xX not because we're using *any* of the
  # xX values but just to initialize xfitX to the desired size (which is the same as the size of xX)
  
  # for 1 to number of timepoints
  for (j in 1:nt) {
    # the expected value of xc(t) at this timepoint t, calculated as a column for all constructs c
    xfit1[,j]<-ac1+fc*time[j]+lambda1[j]
  }
  
  # same as above but for replicate 2
  xfit2<-x2 #for size
  for (j in 1:nt) {
    xfit2[,j]<-ac2+fc*time[j]+lambda2[j]
  }
  
  
  sdfc<-rep(0.1,nx) #standard error of fc
  tstat<-rep(0,nx) # presumably t statistic
  df<-rep(0,nx) # df = vector w one entry for each construct, with value from 0 to -2; 0 if both replicates of relevant construct are good, -1 if just one is, -2 if neither are.  Inits to 0 for all
  p_t<-rep(1,nx) # p value from t test ... initialize to 1 (not significant) for everything
  
  # for 1 to number of constructs
  for (i in 1:nx) {
    # if this construct doesn't have enough "good" measurements, skip it
    if (allbad[i]) next
    
    
    g1<-good1[i,] # g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the first replicate 
    g2<-good2[i,]
    
    # df = vector w one entry for each construct, with value from 0 to -2; 0 if both replicates of relevant construct are good, -1 if just one is, -2 if neither are
    # I suspect that "df" is "degrees of freedom" for each construct
    df[i]<-sum(g1)+sum(g2)-2
    
    # sqrtsum<-function(y) sqrt(sum(y^2))
    # outside of this function, this sdfc value is used only in the output of the construct file
    # I think that sdfc is "standard deviation of fc" for each construct.
    
    # So: while I don't claim to understand why, I know that Roman says:
    # std err of fc = sqrt of sum over t of (lower-case-epsilon for construct c, as a function of t)^2
    # divided by sqrt of (nc - 2)*sum over t of (t^2 - (mean of t)^2)
    # Note that lower-case-epsilon is Xc(t) - xc(t) = xX - xfitX
    # so xfitX - xX = - lower-case-epsilon ... but since it is being squared and then square-rooted, I
    # suppose the negation doesn't matter.
    # nc - 2 is the number of degrees of freedom
    # where nc = number of data points = 2*nt minus any number of points below the threshold
    # (note the description above seems to assume 2 replicates) ... seems to me this must mean 
    # number of data points *for this construct c*, not total.
    # nt in above is number of timepoints, as here ...
    # Roman also says that tc [i.e., the t statistic for construct c] = fc /SE(fc)
    # and the internet tells me that SE(x) = SD(x)/sqrt(n) ... but maybe the sqrt(n) is just the most usual
    # sqrt of degrees of freedom, and could be something else in a more complex system.
    # So, the value being calculated directly below is SD(x), which is why it doesn't have the
    # nc - 2 term that is in the denominator of the SE(x) calculation (in the manuscript, Roman says that
    # fc's sd = sqrt(nc-2)*SE(fc), where nc-2 = degrees of freedom.  Farther below, after the calculation
    # of sdfc, we get the calculation of tc, which is fc/(sdfc/sqrt(df)), where the sdfc/sqrt(df) term 
    # is the calculation of SE(fc).
    
    # result is vector with one std dev of fc for each construct
    sdfc[i]<-sqrtsum( c(xfit1[i,g1],xfit2[i,g2]) - c(x1[i,g1],x2[i,g2]) ) /
      sqrtsum( c(time[g1],time[g2]) - mean(c(time[g1],time[g2])) )
    
  }
  #find median sd
  has_sd<-df>0
  median_sd<-median(sdfc[has_sd])
  sdfc[!has_sd]<-median_sd #just so it isn't 0
  
  # for 1 to number of constructs
  for (i in 1:nx) {
    # don't try to calculate t statistic and p value for any fc that doesn't have a stderr
    if (!has_sd[i]) next
    # calc t statistic of fc
    tstat[i]<-fc[i]/(sdfc[i]/sqrt(df[i]))
    # pt is R function, presumably to do t-test :)  Result must be one-tailed, hence *2
    p_t[i]<-2*pt(-abs(tstat[i]),df=df[i]) #raw p-values from t-test
  }
  
  # ok, lfdr stands for "local fdr", and the lfdr function comes from the qvalue bioconductor package
  # that is installed way above.  The first input is the vector of p-values, and the second input is the
  # estimated proportion of true null p-values, where pi0.method is "the method for automatically choosing tuning
  # parameter in the estimation of π0, the proportion of true null hypotheses."
  lfdr_fc<-rep(1,nx) # nx = number of constructs; default value of lfdr is set to one for all of them
  l<-lfdr(p_t[has_sd],pi0.method="bootstrap")
  lfdr_fc[has_sd]<-l # for constructs that have an sd, the lfdr of the fc could be calculated and is now set
  # to its calculated vale instead of the default
  
  # ac1 is the initial condition (in log2 frequency) for each construct c for replicate 1
  # ac2 is the initial condition (in log2 frequency) for each construct c for replicate 2
  # fc is the fitness of each construct c (calculated across both replicates)
  # sdfc is the std deviation of the fitness of each construct c (calculated across both replicates)
  # p_t is the raw p value of the fc of each construct c (calculated across both replicates)
  # lfdr_fc is the local FDR of each construct (calculated across both replicates)
  # df is the degrees of freedom of each construct c (calculated across both replicates)
  # allbad is a boolean value for each construct c that is true for all the constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments
  vl<-list(ac1,ac2,fc,sdfc,p_t,lfdr_fc,df,allbad)
  return(vl)
}

#iterative robust least squares
# fc is a symmetric probe-by-probe matrix where the value in each cell is the fc value for the construct including
# those two probes, or zero if there is no construct for that probe-by-probe combination.
# w0 is a symmetric probe-by-probe matrix of initial weights, where the value in each cell is zero for "bad"
# and/or non-existent constructs and one for "good" constructs.
# The manuscript gives eqn 2 as:
# fc = fp + fp' + pi sub pp'
# and states "The probe fitnesses are found by robust fitting of Eq. (2). 
# The probe-level pi-scores are the residuals of the robust fit."
irls<-function(fc,w0,probes,ag=2,tol=1e-3,maxiter=30) {
  # w0 is the physical goodness of constructs. It is not subject to change. # ab: what is "physical goodness"??
  # It is used to modify w, to silence bad constructs 
  
  # upper.tri seems to be getting the upper triangle of the symmetric probe-by-probe fc matrix
  # What data structure is expressed_utri? a matrix of fc values?  or a matrix of booleans? or ...?
  expressed_utri<-upper.tri(fc) & w0>0
  n<-dim(fc)[1] # lessee--that would be number of probes; note that n here is *different* than n outside the scope of this function, where it is the number of genes :(
  w<-matrix(1,nrow=n,ncol=n) #initial weights # ab: probe-by-probe matrix with default values of 1 everywhere except on diagonal, where they are zero
  diag(w)<-0
  # errr.  what are e and f?
  fij<-matrix(0,nrow=n,ncol=n) #initial weights
  eij<-matrix(0,nrow=n,ncol=n) #initial weights
  b<-rep(0,n) #rhs #ab: wth is rhs?
  
  #iteration step 0
  w<-w*w0 # I think here we are "silenc[ing] bad constructs", which have w0 = 0
  A<-w
  for (i in 1:n) { # for each probe
    b[i]<-sum(fc[,i]*w[,i]) # sum of fc*weight for all probes paired with this probe.  Why set i as the second index into in the probe-by-probe matrix rather than the first?  Isn't the matrix symmetric?
    A[i,i]<-sum(w[,i])+small # reminder: small<-1e-6 ; it is a parameter hardcoded at the top of the whole notebook.
    # so, A = w and w is 1 everywhere except on the diagonal and at probe-by-probe pairs where the analogous 
    # construct is either bad or non-existent.  This seems to be setting the values along the diagonal to the
    # sum of the weights for that probe (across all other probes it is paired with) plus a tiny value, presumably
    # to make sure the diagonal is now *never* zero.
  }
  
  # apparently, "solve() function solves equation a %*% x = b for x, where b is a vector or matrix." (from http://www.endmemo.com/program/R/solve.php)
  # Note that %*% is apparently R notation for "multiply matrices".  Ok, so figure out what matrix you need to multiply A by to get b,
  # and call that matrix y (again, note y here is *different* from y outside this function, where it is #log2 frequencies :( )
  # Here, y is a matrix with 1 column and as many rows as there are probes, since b is an array (essentially, a matrix with 1 row and as many columns as there are probes)
  y<-solve(A,b)
  names(y)<-probes
  # ok, here's that same logic for getting all pairs of probes
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      fij[i,j]<-y[i]+y[j]
    }
  }
  fij<-fij+t(fij) # make fij symmetric
  eij<-fc-fij #residuals
  # ab: end iteration step 0
  
  logStrings = vector('character') # starts empty
  l<-1 #counter
  rel<-1 # apparently, to judge by Roman's comment at the end of this while loop, "rel" = "relative error"
  while (rel > tol & l < maxiter) {#iterate until tol is reached or maxit
    s<-sd(eij[expressed_utri]) #calculate sd only from expressed constructs # ab: calc stddev of residuals for all expressed probes
    yold<-y # archive y, the matrix w/ one column and rows for each probe
    
    # I see that ag = 2 (as passed in to function params) but wth *is* ag?
    # and ... just *what* is this w calculation?  Where does it come from?  What does it *mean*?
    w<-(1-eij^2/(ag*s)^2)^2 #something like Tukey's biweight
    w[abs(eij)>(ag*s)]<-0
    diag(w)<-0
    
    #Eeep! this is an exact repeat of the code above in the section labeled "iteration step 0", except that the
    # fij<-matrix(0,nrow=n,ncol=n) has been moved inside the iteration instead of being before it.
    w<-w*w0
    
    A<-w
    for (i in 1:n) {
      b[i]<-sum(fc[,i]*w[,i])
      A[i,i]<-sum(w[,i])+small
    }
    y<-solve(A,b)
    names(y)<-probes
    fij<-matrix(0,nrow=n,ncol=n) #initial weights
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        fij[i,j]<-y[i]+y[j]
      }
    }
    fij<-fij+t(fij)
    eij<-fc-fij #residuals
    
    # ok, in iteration 0, rel = one and l (the counter) = one
    # here, the counter = the counter plus one
    # rel: we're calculating the difference of the y value for this time through the loop
    # from the y value the last time through the loop and then squaring it (presumably to get rid of 
    # any negatives); we do this across all probes and sum the squared values.  Then we divide that sum
    # of squares by either one or the sum across all probes of the square of the y values from the last time 
    # through the loop, whichever is greater.  Then we take the sqrt of that ratio and say it is the 
    # "relative error" ... whatever that means.
    
    rel<-sqrt(sum((y-yold)^2)/max(1,sum(yold^2))) #relative error
    # reminder: sqrtsum<-function(y) sqrt(sum(y^2))
    logStrings= append(logStrings, capture.output(cat(l,sqrtsum(yold),sqrtsum(y-yold),"\n"))) # print out some values to the screen
    l<-l+1
    
    # we break out of this loop when either: 
    # a) the relative error has fallen below the "tol" value (what does "tol" stand for?) OR
    # b) we have performed the maximum number of trips through the loop.
    # Both tol and maxiter are parameters passed into this function.
  }
  # ok, would it kill ya to give me a hint of what these per-probe values in y *are*?
  # Roman's comment when y is unpacked from this returned list is "#these are probe fitnesses fp"
  # Roman's comment when eij is unpacked from this returned list is "#raw pi-scores per construct"
  vl<-list(y, fij, eij, logStrings)
  return(vl)
}

calcFdrsLeftAndRight<-function(pi_mean, pi_null){
  enull<-ecdf(pi_null) # used in plot of pi scores by fdr
  emean<-ecdf(pi_mean) # used in plot of pi scores by fdr  
  
  fdr_left<-pmin(1,enull(pi_mean)/emean(pi_mean))
  fdr_right<-pmin(1,(enull(-pi_mean))/(1-emean(pi_mean)))
  return(list(fdr_left=fdr_left, fdr_right=fdr_right))
}

# TODO: are these bads for the 2 genes or the 2 replicates?  Former ok, latter not
# so ... this value is assigned here but then never used in anything that comes after it; commenting it out to see if anything breaks :)
# good<-!bad1 & !bad2

prepData=function(input_filename, nt){
  X<-read.table(input_filename,sep="\t",header=TRUE)
  #preliminary preparations of the input data frame
  # TODO: ok, the 5 & 6 are definitely hardcoding for #s of expected info columns--perhaps not a problem.
  # TODO: the 2 seems like it is probably hardcoding for # of replicates, and thus a problem
  data<-data.matrix(X[,6:(5+2*nt)])
  # TODO: geneA and geneB are hardcodings for info column names; undesirable.
  # TODO: however, assuming that there will be two genes (and two probes) per construct is acceptable
  # since this is a DUAL crispr pipeline.
  good<-(X$geneA != X$geneB) #reject any constructs with two 0's
  goodX<-X[good,] #the 0-0 constructs are gone
  nn<-sum(good) #this many constructs
  
  cpA<-as.character(goodX$probeA)
  # TODO: "Nontargeting" (throughout) is hardcoding for ntc prefix; should be refactored
  ix<-grep("NonTargeting",cpA)
  cpA[ix]<-paste("0",cpA[ix],sep="") #this puts NonTargeting probes at the beginning of alphabetically sorted order
  
  cpB<-as.character(goodX$probeB)
  ix<-grep("NonTargeting",cpB)
  cpB[ix]<-paste("0",cpB[ix],sep="")
  
  pswitch<-cpA>cpB #need to switch?
  phold<-cpA[pswitch]
  cpA[pswitch]<-cpB[pswitch]
  cpB[pswitch]<-phold #cpA and cpB are always in alphabetical order, cpA < cpB
  probes<-sort(unique(c(cpA,cpB))) #entire probe set in alphabetical order
  nprobes<-length(probes)
  
  
  cgA<-as.character(goodX$geneA)
  cgB<-as.character(goodX$geneB)
  genes<-sort(unique(cgA)) #should be 74 "genes"
  n<-length(genes) # n = 74 if doing it by genes or 222 if doing it by probe # number of genes
  mm<-n*(n-1)/2
  
  gswitch<-cgA>cgB #need to switch?
  ghold<-cgA[gswitch]
  cgA[gswitch]<-cgB[gswitch]
  cgB[gswitch]<-ghold
  
  # TODO: probably separator should be specified in set-up rather than hardcoded
  gA_gB<-paste(cgA,cgB,sep="_")
  pA_pB<-paste(cpA,cpB,sep="_")
  goodX<-data.frame(goodX,cgA,cgB,gA_gB) #now gA_gB is ordered so that gA < gB
  
  # TODO: def refactor here, as above
  gooddata<-data.matrix(goodX[,6:(5+2*nt)])
  # TODO: are we ok with hardcoding what the pseudocount will be?  Maybe should set in params.
  gooddata[gooddata==0]<-1 #pseudocounts
  abundance<-apply(gooddata,2,sum)
  y<-t(log2(t(gooddata)/abundance)) #log2 frequencies
  #end data prep 
  
  return(list(nn=nn, probes=probes, nprobes=nprobes, genes=genes, n=n, pA_pB=pA_pB, y=y))  
}

getXsAbsAndBads=function(y, nt, ab0) {
  # for replicate 1
  # "2"s here and line below look like replicate hardcoding; refactor
  
  # Ok, so this is getting the log2 frequencies for the 1st replicate of all timepts
  x1<-y[,seq(1,2*nt,by=2)]
  # This is getting the abundance thresholds for all of these 1st replicates
  ab1<-ab0[seq(1,2*nt,by=2)]
  # TODO: 1st "2" here is ok (designates row vs col) ... but what is source of second?
  bad1<-apply(t(x1)>ab1,2,sum)<2
  
  # for replicate 2
  # "2"s here and line below look like replicate hardcoding; refactor
  # Ok, so this is getting the log2 frequencies for the 2nd replicate of all timepts
  x2<-y[,seq(2,2*nt,by=2)]
  # This is getting the abundance thresholds for all of these 2nd replicates
  ab2<-ab0[seq(2,2*nt,by=2)]
  # TODO: 1st "2" here is ok (designates row vs col); 2nd doesn't mean replicates but 
  # represents arbitrary(?) decision that a construct is bad if it isn't over the relevant
  # abundance threshold in at least two timepoints.  Could one reasonably choose a 
  # different value?  If so, should refactor to set in params.
  
  # 2 in apply call params means "sum over columns".  Because of t call, x2 has been 
  # transposed so t(x2) has timept as rows and constructs as columns.
  # We're calling "bad" any construct that doesn't have log2 frequency values above the timepoint-specific
  # abundance threshold for at least two timepoints
  bad2<-apply(t(x2)>ab2,2,sum)<2
  return(list(x1=x1, ab1=ab1, bad1=bad1, x2=x2, ab2=ab2, bad2=bad2))
}
