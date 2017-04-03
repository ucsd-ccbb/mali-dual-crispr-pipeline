#scalar version
Var <- function(x)
  mean(x ^ 2) - mean(x) ^ 2

Cov <- function(x, y)
  mean(x * y) - mean(x) * mean(y)

sqrtsum <- function(y)
  sqrt(sum(y ^ 2))

temp4 <- function(i, timepointsList, good1, x1, good2, x2){
  ac1 <- x1[i, 1]
  ac2 <- x2[i, 1] #just a guess
  fc <- 0 # nx = number of constructs  
  
  # apparently "f" stands for fitness--which is covariance (see below)--and # stands for replicate
  # v stands for variance and # " "
  f1 <- 0
  f2 <- 0
  v1 <- 0
  v2 <- 0
  
  # g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the first replicate
  g1 <- good1[i,]     
  if (sum(g1) > 1) {
    #it's a good experiment # if there are at least two good timepoints for this construct in replicate 1
    
    # get the mean of the log2 frequencies for all the good timepoints for this construct in replicate 1
    # timepointsList is GLOBAL vector of timepoints (numbers, usually in days, in ascending order)
    mx1 <- mean(x1[i, g1]) 
    # get mean of timepoints (e.g. mean number of days) for all the good timepoints for this construct in replicate 1
    mt1 <- mean(timepointsList[g1])       
    v1 <- Var(timepointsList[g1]) # get variance of timepoints " "
    # f1 = covariance of log2 frequencies for all the good timepoints for this construct in replicate 1 with the timepoints for those good timepoints
    f1 <- Cov(x1[i, g1], timepointsList[g1])
  }
  
  # do the exact same thing as above, but for replicate 2
  g2 <- good2[i,]
  if (sum(g2) > 1) {
    #it's a good experiment
    mx2 <- mean(x2[i, g2])
    mt2 <- mean(timepointsList[g2])
    v2 <- Var(timepointsList[g2])
    f2 <- Cov(x2[i, g2], timepointsList[g2])
  }
  
  # fc[i] is the combined fitness (across replicates) for construct i .
  #the combined fitness from replicate 1+2
  #fc remains defined up to an additive constant
  fc <- (f1 + f2) / (v1 + v2) 
  
  if (sum(g1) > 1) {
    # if there are at least two good timepoints for this construct in replicate 1
    # ac1 = mean of log2 freqs for this construct for good timepoints for rep 1 - (mean of timepts for good timepoints for this construct for rep 1)*combined fitness across replicates for this construct
    ac1 <- mx1 - fc * mt1
  }
  
  # same as above but for replicate 2
  if (sum(g2) > 1) {
    ac2 <- mx2 - fc * mt2
  }  
  
  return(list(fc = fc, ac1 = ac1, ac2 = ac2))
}

temp1 <- function(x1, ab1, x2, ab2, nt, nx, timepointsList){
  # Prediction: useless1 and useless2 here should have the same values as bad1 and bad2 in the
  # Number of Constructs Below Abundance Threshold section ... boxed code seems to calc same result
  # as code in that section, but by slightly diff (but equivalent) means
  # -------------------
  # t(x1) has timept as rows and constructs as columns.
  # I can't for the life of me figure out why it is necessary to transpose x1, compare it to ab1,
  # and then transpose it back ... since the ">ab1" evaluation will be at the individual cell level,
  # and the individual cells hold log2 frequency values no matter which way you set up the rows and cols,
  # and since no margin calculations are being done in the good1 calculation, .... WHY?
  
  # constructs as rows, timept for 1st replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2 freq for row/col combination is above relevant abundance threshold
  good1 <- t(t(x1) > ab1) 
  # constructs as rows, timept for 2nd replicates ONLY as cols, cell values are 0/FALSE, 1/TRUE for whether log2 freq for row/col combination is above relevant abundance threshold
  good2 <- t(t(x2) > ab2)   
  # 1 = sum over rows--i.e., constructs.  Any construct that isn't above abundance threshold in 1st replicate in at least 2 timepoints has 1/TRUE in "useless" matrix
  useless1 <- apply(good1, 1, sum) < 2     
  useless2 <- apply(good2, 1, sum) < 2
  
  # ab
  #print(sum(useless1))
  #print(sum(useless2))
  # PREDICTION CORRECT!
  # --------------------
  
  # I do not understand how the statements below actually *remove* anything from the good1 and good2 dfs ...
  # It seems like, for each constructs that isn't above abundance thresholds in at least two timepoints for this replicate, it sets that construct's "good" values to FALSE for *all* timepoints in this replicate
  good1[useless1,] <- FALSE #remove singletons
  good2[useless2,] <- FALSE #remove singletons
  
  # Note: here apply(goodX,1,sum) is NOT uselessX because goodX was changed above
  #allbad is true for all the constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments
  #in this case I have nothing to use in either experiment
  allbad <- apply(good1, 1, sum) < 2 & apply(good2, 1, sum) < 2 
  
  #just a guess # log2 frequencies for all constructs for first timepoint in this replicate
  ac1 <- x1[, 1]
  ac2 <- x2[, 1] #just a guess
  fc <- rep(0, nx) # nx = number of constructs
  
  # for 1 to number of constructs
  for (i in 1:nx) {
    # if this construct doesn't have at least two timepoints above abundance threshold in at least one replicate, ignore it and move on
    if (allbad[i])
      next #from now on there is at least one good experiment
    
    temp4Results = temp4(i, timepointsList, good1, x1, good2, x2)
    ac1[i] = temp4Results$ac1
    ac2[i] = temp4Results$ac2
    fc[i] = temp4Results$fc
  }  
  
  # ac is the initial condition (in log2 frequency) for construct c
  # this is normalizing ac1 ... I think this is what is going on:
  # Roman's methods say "By definition, log2 relative frequencies satisfy the constraint
  # sum over c of (2^xc) = 1 at all times"
  # ac1 is the set of log2 frequencies for all constructs in replicate 1 at "initial conditions", so
  # ac1 must satisfy the above constraint.  IFF the above constraint is satisfied, -log2(1) = 0.
  # If the above constraint isn't satisfied, the value of -log2(sum over c of (2^ac)) is subtracted from
  # ac to ensure the constraint is satisfied.
  alpha <- -log2(sum(2 ^ ac1))
  ac1 <- ac1 + alpha #enforce normalization at time=0, sum(2^ac)=1
  
  # same as above but for replicate 2
  alpha <- -log2(sum(2 ^ ac2))
  ac2 <- ac2 + alpha #enforce normalization at time=0, sum(2^ac)=1  
  
  return(list(ac1 = ac1, ac2 = ac2, fc = fc, allbad = allbad, good1 = good1, good2 = good2))
}

temp2 <- function(nx, has_sd, fc, df, sdfc) {
  library(qvalue)
  
  l <- 0
  tstat <- rep(0, nx) # presumably t statistic
  # p value from t test ... initialize to 1 (not significant) for everything
  p_t <- rep(1, nx)   
  
  # for 1 to number of constructs
  for (i in 1:nx) {
    # don't try to calculate t statistic and p value for any fc that doesn't have a stderr
    if (!has_sd[i])
      next
    # calc t statistic of fc
    tstat[i] <- fc[i] / (sdfc[i] / sqrt(df[i]))
    # pt is R function, presumably to do t-test :)  Result must be one-tailed, hence *2
    #raw p-values from t-test
    p_t[i] <- 2 * pt(-abs(tstat[i]), df = df[i]) 
  }
  
  # ok, lfdr stands for "local fdr", and the lfdr function comes from the qvalue bioconductor package
  # that is installed way above.  The first input is the vector of p-values, and the second input is the
  # estimated proportion of true null p-values, where pi0.method is "the method for automatically choosing tuning
  # parameter in the estimation of π0, the proportion of true null hypotheses."
  
  # nx = number of constructs; default value of lfdr is set to one for all of them
  lfdr_fc <- rep(1, nx) 
  l <- lfdr(p_t[has_sd], pi0.method = "bootstrap")
  
  # for constructs that have an sd, the lfdr of the fc could be calculated and is now set
  # to its calculated vale instead of the default    
  lfdr_fc[has_sd] <- l 
  
  return(list(p_t = p_t, lfdr_fc = lfdr_fc))
}

temp3 <- function(nt, timepointsList, ac1, ac2, fc, x1, x2, nx, allbad, good1, good2) {
  # ok, the *expected* log2 frequency for construct x at time t is:
  # xc(t) = ac + fc*t - log2(sum over c of 2^(ac + fc*t))
  # apparently the below code calls the term starting with "-log2" "lambda".
  # Note that lambda is being calculated separately for each combination of timepoint+replicate
  
  # nt = number of timepoints
  lambda1 <- rep(0, nt)
  lambda2 <- rep(0, nt)
  
  # for 1 to number of timepoints
  for (i in 1:nt) {
    lambda1[i] <- -log2(sum(2 ^ (ac1 + fc * timepointsList[i])))
    lambda2[i] <- -log2(sum(2 ^ (ac2 + fc * timepointsList[i])))
  } #these are initial estimates of lambda(t)
  
  # x1 is log2 frequencies for the 1st replicate of all timepts
  #for size # I think this means that xfitX is being set to xX not because we're using *any* of the
  # xX values but just to initialize xfitX to the desired size (which is the same as the size of xX)
  xfit1 <- x1
  
  # for 1 to number of timepoints
  for (j in 1:nt) {
    # the expected value of xc(t) at this timepoint t, calculated as a column for all constructs c
    xfit1[, j] <- ac1 + fc * timepointsList[j] + lambda1[j]
  }
  
  # same as above but for replicate 2
  xfit2 <- x2 #for size
  for (j in 1:nt) {
    xfit2[, j] <- ac2 + fc * timepointsList[j] + lambda2[j]
  }
  
  sdfc <- rep(0.1, nx) #standard error of fc
  # df = vector w one entry for each construct, with value from 0 to -2; 0 if both replicates of relevant construct are good, -1 if just one is, -2 if neither are.  Inits to 0 for all
  df <- rep(0, nx)     
  
  # for 1 to number of constructs
  for (i in 1:nx) {
    # if this construct doesn't have enough "good" measurements, skip it
    if (allbad[i])
      next
    
    # g1 = true/false values of whether construct i passes various abundance filters for all timepoints for the first replicate
    g1 <- good1[i,]       
    g2 <- good2[i,]
    
    # df = vector w one entry for each construct, with value from 0 to -2; 0 if both replicates of relevant construct are good, -1 if just one is, -2 if neither are
    # I suspect that "df" is "degrees of freedom" for each construct
    df[i] <- sum(g1) + sum(g2) - 2
    
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
    sdfc[i] <-
      sqrtsum(c(xfit1[i, g1], xfit2[i, g2]) - c(x1[i, g1], x2[i, g2])) /
      sqrtsum(c(timepointsList[g1], timepointsList[g2]) - mean(c(timepointsList[g1], timepointsList[g2])))
    
  }
  #find median sd
  has_sd <- df > 0
  median_sd <- median(sdfc[has_sd])
  sdfc[!has_sd] <- median_sd #just so it isn't 0
  
  return(list(df = df, has_sd = has_sd, sdfc = sdfc))
}

# x1 is log2 frequencies for the 1st replicate of all timepts
# x2 is log2 frequencies for the 2nd replicate of all timepts
# ab1 is abundance thresholds for all 1st replicates
# ab2 is abundance thresholds for all 2nd replicates
fit_ac_fc <- function(x1, ab1, x2, ab2, nt, timepointsList) {
  nx <- nrow(x1)
  
  tempResult = temp1(x1, ab1, x2, ab2, nt, nx, timepointsList)
  ac1 = tempResult$ac1
  ac2 = tempResult$ac2
  fc = tempResult$fc
  allbad = tempResult$allbad
  good1 = tempResult$good1
  good2 = tempResult$good2

  temp3Results = temp3(nt, timepointsList, ac1, ac2, fc, x1, x2, nx, allbad, good1, good2)
  df = temp3Results$df
  has_sd = temp3Results$has_sd
  sdfc = temp3Results$sdfc
  
  temp2Results = temp2(nx, has_sd, fc, df, sdfc) 
  p_t = temp2Results$p_t
  lfdr_fc = temp2Results$lfdr_fc
  
  # ac1 is the initial condition (in log2 frequency) for each construct c for replicate 1
  # ac2 is the initial condition (in log2 frequency) for each construct c for replicate 2
  # fc is the fitness of each construct c (calculated across both replicates)
  # sdfc is the std deviation of the fitness of each construct c (calculated across both replicates)
  # p_t is the raw p value of the fc of each construct c (calculated across both replicates)
  # lfdr_fc is the local FDR of each construct (calculated across both replicates)
  # df is the degrees of freedom of each construct c (calculated across both replicates)
  # allbad is a boolean value for each construct c that is true for all the constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments
  vl <- list(ac1, ac2, fc, sdfc, p_t, lfdr_fc, df, allbad)
  return(vl)
}

calcFdrsLeftAndRight <- function(pi_mean, pi_null) {
  enull <- ecdf(pi_null) # used in plot of pi scores by fdr
  emean <- ecdf(pi_mean) # used in plot of pi scores by fdr
  
  fdr_left <- pmin(1, enull(pi_mean) / emean(pi_mean))
  fdr_right <- pmin(1, (enull(-pi_mean)) / (1 - emean(pi_mean)))
  return(list(fdr_left = fdr_left, fdr_right = fdr_right))
}

prepData = function(input_filename, nt) {
  X <- read.table(input_filename, sep = "\t", header = TRUE)
  #preliminary preparations of the input data frame
  # TODO: ok, the 5 & 6 are definitely hardcoding for #s of expected info columns--perhaps not a problem.
  # TODO: the 2 seems like it is probably hardcoding for # of replicates, and thus a problem
  #data <- data.matrix(X[, 6:(5 + 2 * nt)])
  # TODO: geneA and geneB are hardcodings for info column names; undesirable.
  # TODO: however, assuming that there will be two genes (and two probes) per construct is acceptable
  # since this is a DUAL crispr pipeline.
  good <- (X$geneA != X$geneB) #reject any constructs with two 0's
  goodX <- X[good,] #the 0-0 constructs are gone
  nn <- sum(good) #this many constructs
  
  cpA <- as.character(goodX$probeA)
  # TODO: "Nontargeting" (throughout) is hardcoding for ntc prefix; should be refactored
  ix <- grep("NonTargeting", cpA)
  cpA[ix] <-
    paste("0", cpA[ix], sep = "") #this puts NonTargeting probes at the beginning of alphabetically sorted order
  
  cpB <- as.character(goodX$probeB)
  ix <- grep("NonTargeting", cpB)
  cpB[ix] <- paste("0", cpB[ix], sep = "")
  
  pswitch <- cpA > cpB #need to switch?
  phold <- cpA[pswitch]
  cpA[pswitch] <- cpB[pswitch]
  cpB[pswitch] <-
    phold #cpA and cpB are always in alphabetical order, cpA < cpB
  probes <-
    sort(unique(c(cpA, cpB))) #entire probe set in alphabetical order
  nprobes <- length(probes)
  
  
  cgA <- as.character(goodX$geneA)
  cgB <- as.character(goodX$geneB)
  genes <- sort(unique(cgA)) #should be 74 "genes"
  n <-
    length(genes) # n = 74 if doing it by genes or 222 if doing it by probe # number of genes
  #mm <- n * (n - 1) / 2
  
  gswitch <- cgA > cgB #need to switch?
  ghold <- cgA[gswitch]
  cgA[gswitch] <- cgB[gswitch]
  cgB[gswitch] <- ghold
  
  # TODO: probably separator should be specified in set-up rather than hardcoded
  gA_gB <- paste(cgA, cgB, sep = "_")
  pA_pB <- paste(cpA, cpB, sep = "_")
  goodX <-
    data.frame(goodX, cgA, cgB, gA_gB) #now gA_gB is ordered so that gA < gB
  
  # TODO: def refactor here, as above
  gooddata <- data.matrix(goodX[, 6:(5 + 2 * nt)])
  # TODO: are we ok with hardcoding what the pseudocount will be?  Maybe should set in params.
  gooddata[gooddata == 0] <- 1 #pseudocounts
  abundance <- apply(gooddata, 2, sum)
  y <- t(log2(t(gooddata) / abundance)) #log2 frequencies
  #end data prep
  
  return(list(
    nn = nn,
    probes = probes,
    nprobes = nprobes,
    genes = genes,
    n = n,
    pA_pB = pA_pB,
    y = y
  ))
}

getXsAbsAndBads = function(y, nt, ab0) {
  # for replicate 1
  # "2"s here and line below look like replicate hardcoding; refactor
  # Ok, so this is getting the log2 frequencies for the 1st replicate of all timepts
  x1 <- y[, seq(1, 2 * nt, by = 2)]
  # This is getting the abundance thresholds for all of these 1st replicates
  ab1 <- ab0[seq(1, 2 * nt, by = 2)]
  # TODO: 1st "2" here is ok (designates row vs col) ... but what is source of second?
  bad1 <- apply(t(x1) > ab1, 2, sum) < 2
  
  # for replicate 2
  # "2"s here and line below look like replicate hardcoding; refactor
  # Ok, so this is getting the log2 frequencies for the 2nd replicate of all timepts
  x2 <- y[, seq(2, 2 * nt, by = 2)]
  # This is getting the abundance thresholds for all of these 2nd replicates
  ab2 <- ab0[seq(2, 2 * nt, by = 2)]
  # TODO: 1st "2" here is ok (designates row vs col); 2nd doesn't mean replicates but
  # represents arbitrary(?) decision that a construct is bad if it isn't over the relevant
  # abundance threshold in at least two timepoints.  Could one reasonably choose a
  # different value?  If so, should refactor to set in params.
  
  # 2 in apply call params means "sum over columns".  Because of t call, x2 has been
  # transposed so t(x2) has timept as rows and constructs as columns.
  # We're calling "bad" any construct that doesn't have log2 frequency values above the timepoint-specific
  # abundance threshold for at least two timepoints
  bad2 <- apply(t(x2) > ab2, 2, sum) < 2
  return(list(
    x1 = x1,
    ab1 = ab1,
    bad1 = bad1,
    x2 = x2,
    ab2 = ab2,
    bad2 = bad2
  ))
}
