small <- 1e-6

#iterative robust least squares
# fc is a symmetric probe-by-probe matrix where the value in each cell is the fc value for the construct including
# those two probes, or zero if there is no construct for that probe-by-probe combination.
# w0 is a symmetric probe-by-probe matrix of initial weights, where the value in each cell is zero for "bad"
# and/or non-existent constructs and one for "good" constructs.
# The manuscript gives eqn 2 as:
# fc = fp + fp' + pi sub pp'
# and states "The probe fitnesses are found by robust fitting of Eq. (2).
# The probe-level pi-scores are the residuals of the robust fit."
irls <- function(fc,
                 w0,
                 probes,
                 ag = 2,
                 tol = 1e-3,
                 maxiter = 30) {
  # w0 is the physical goodness of constructs. It is not subject to change. # ab: what is "physical goodness"??
  # It is used to modify w, to silence bad constructs
  
  # upper.tri seems to be getting the upper triangle of the symmetric probe-by-probe fc matrix
  # What data structure is expressed_utri? a matrix of fc values?  or a matrix of booleans? or ...?
  expressed_utri <- upper.tri(fc) & w0 > 0
  n <-
    dim(fc)[1] # lessee--that would be number of probes; note that n here is *different* than n outside the scope of this function, where it is the number of genes :(
  w <-
    matrix(1, nrow = n, ncol = n) #initial weights # ab: probe-by-probe matrix with default values of 1 everywhere except on diagonal, where they are zero
  diag(w) <- 0
  # errr.  what are e and f?
  fij <- matrix(0, nrow = n, ncol = n) #initial weights
  eij <- matrix(0, nrow = n, ncol = n) #initial weights
  b <- rep(0, n) #rhs #ab: wth is rhs?
  
  #iteration step 0
  w <-
    w * w0 # I think here we are "silenc[ing] bad constructs", which have w0 = 0
  A <- w
  for (i in 1:n) {
    # for each probe
    b[i] <-
      sum(fc[, i] * w[, i]) # sum of fc*weight for all probes paired with this probe.  Why set i as the second index into in the probe-by-probe matrix rather than the first?  Isn't the matrix symmetric?
    A[i, i] <-
      sum(w[, i]) + small # reminder: small<-1e-6 ; it is a parameter hardcoded at the top of the whole notebook.
    # so, A = w and w is 1 everywhere except on the diagonal and at probe-by-probe pairs where the analogous
    # construct is either bad or non-existent.  This seems to be setting the values along the diagonal to the
    # sum of the weights for that probe (across all other probes it is paired with) plus a tiny value, presumably
    # to make sure the diagonal is now *never* zero.
  }
  
  # apparently, "solve() function solves equation a %*% x = b for x, where b is a vector or matrix." (from http://www.endmemo.com/program/R/solve.php)
  # Note that %*% is apparently R notation for "multiply matrices".  Ok, so figure out what matrix you need to multiply A by to get b,
  # and call that matrix y (again, note y here is *different* from y outside this function, where it is #log2 frequencies :( )
  # Here, y is a matrix with 1 column and as many rows as there are probes, since b is an array (essentially, a matrix with 1 row and as many columns as there are probes)
  y <- solve(A, b)
  names(y) <- probes
  # ok, here's that same logic for getting all pairs of probes
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      fij[i, j] <- y[i] + y[j]
    }
  }
  fij <- fij + t(fij) # make fij symmetric
  eij <- fc - fij #residuals
  # ab: end iteration step 0
  
  logStrings = vector('character') # starts empty
  l <- 1 #counter
  rel <-
    1 # apparently, to judge by Roman's comment at the end of this while loop, "rel" = "relative error"
  while (rel > tol &
         l < maxiter) {
    #iterate until tol is reached or maxit
    s <-
      sd(eij[expressed_utri]) #calculate sd only from expressed constructs # ab: calc stddev of residuals for all expressed probes
    yold <-
      y # archive y, the matrix w/ one column and rows for each probe
    
    # I see that ag = 2 (as passed in to function params) but wth *is* ag?
    # and ... just *what* is this w calculation?  Where does it come from?  What does it *mean*?
    w <- (1 - eij ^ 2 / (ag * s) ^ 2) ^ 2 #something like Tukey's biweight
    w[abs(eij) > (ag * s)] <- 0
    diag(w) <- 0
    
    #Eeep! this is an exact repeat of the code above in the section labeled "iteration step 0", except that the
    # fij<-matrix(0,nrow=n,ncol=n) has been moved inside the iteration instead of being before it.
    w <- w * w0
    
    A <- w
    for (i in 1:n) {
      b[i] <- sum(fc[, i] * w[, i])
      A[i, i] <- sum(w[, i]) + small
    }
    y <- solve(A, b)
    names(y) <- probes
    fij <- matrix(0, nrow = n, ncol = n) #initial weights
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        fij[i, j] <- y[i] + y[j]
      }
    }
    fij <- fij + t(fij)
    eij <- fc - fij #residuals
    
    # ok, in iteration 0, rel = one and l (the counter) = one
    # here, the counter = the counter plus one
    # rel: we're calculating the difference of the y value for this time through the loop
    # from the y value the last time through the loop and then squaring it (presumably to get rid of
    # any negatives); we do this across all probes and sum the squared values.  Then we divide that sum
    # of squares by either one or the sum across all probes of the square of the y values from the last time
    # through the loop, whichever is greater.  Then we take the sqrt of that ratio and say it is the
    # "relative error" ... whatever that means.
    
    rel <- sqrt(sum((y - yold) ^ 2) / max(1, sum(yold ^ 2))) #relative error
    # reminder: sqrtsum<-function(y) sqrt(sum(y^2))
    logStrings = append(logStrings, capture.output(cat(
      l, sqrtsum(yold), sqrtsum(y - yold), "\n"
    ))) # print out some values to the screen
    l <- l + 1
    
    # we break out of this loop when either:
    # a) the relative error has fallen below the "tol" value (what does "tol" stand for?) OR
    # b) we have performed the maximum number of trips through the loop.
    # Both tol and maxiter are parameters passed into this function.
  }
  # ok, would it kill ya to give me a hint of what these per-probe values in y *are*?
  # Roman's comment when y is unpacked from this returned list is "#these are probe fitnesses fp"
  # Roman's comment when eij is unpacked from this returned list is "#raw pi-scores per construct"
  vl <- list(y, fij, eij, logStrings)
  return(vl)
}

doRobustFitting <-
  function(nn,
           pA_pB,
           allbad,
           fc,
           nprobes,
           probes,
           n,
           pp_fc,
           sdfc, genes, niter) {
    # nn = number of *good* constructs
    # whatever u is, it is 0 for all bad and/or non-existent constructs and, at least at the beginning, one for all good constructs
    u1 <- rep(0, nn)
    names(u1) <- pA_pB
    u1[!allbad] <- 1  #all other weights set to 1
    
    # argh, heaven preserve me from variables fc0 and fc_0 that mean different things :(
    fc0 <- fc
    
    # just sets expected matrix size for all these--num probes by num probes, w/defaults of zero, and row/col names
    fc_0 <- matrix(0, nrow = nprobes, ncol = nprobes)
    sdfc_0 <- matrix(0, nrow = nprobes, ncol = nprobes)
    w0_0 <- matrix(0, nrow = nprobes, ncol = nprobes)
    pp_0 <- matrix(0, nrow = nprobes, ncol = nprobes)
    rownames(fc_0) <- probes
    colnames(fc_0) <- probes
    rownames(sdfc_0) <- probes
    colnames(sdfc_0) <- probes
    rownames(w0_0) <- probes
    colnames(w0_0) <- probes
    rownames(pp_0) <- probes
    colnames(pp_0) <- probes
    
    
    # for 1 to number of genes - 1
    for (i in 1:(n - 1)) {
      # TODO: Ok, where did these "3"s come from?  From the expected number of probes per gene?
      # If so, that's *got* to be refactored out ...
      for (k in 1:3) {
        # index of current probe for current gene in probes array
        iprobe <- (i - 1) * 3 + k
        # probes is entire probe name set in alphabetical order
        iprobe_name <- probes[iprobe]
        
        # note: for (i in 1:(n-1)) {
        #for (j in (i+1):n) {
        # nested loops here are the same as used above to generate ggnames...
        # for all possible second genes that come after the current first gene (i)
        for (j in (i + 1):n) {
          # for all the constructs of the second gene
          for (l in 1:3) {
            jprobe <- (j - 1) * 3 + l
            jprobe_name <-
              probes[jprobe] # as for first gene, get index and then name of second gene
            
            #generate the construct name
            construct <- paste(iprobe_name, "_", jprobe_name, sep = "")
            w0_0[iprobe_name, jprobe_name] <-
              u1[construct] #initial weights. non-existent pairs will have w0=0.  #ab: Roman's comment here says this is initial weights, but his comment in the irls function says that this is actually the "physical goodness" of each construct.  Seems to be integer booleans (i.e., 1 for true, 0 for false)
            pp_0[iprobe_name, jprobe_name] <-
              pp_fc[construct] #initial weights. non-existent pairs will have w0=0 # ab: I don't think this comment is right either.  pp_fc are the posterior probabilities for the fc of each construct
            fc_0[iprobe_name, jprobe_name] <- fc0[construct]
            sdfc_0[iprobe_name, jprobe_name] <- sdfc[construct]
          }
        }
      }
    }
    w0_0 <- w0_0 + t(w0_0) #make symmetric
    fc_0 <- fc_0 + t(fc_0)
    pp_0 <- pp_0 + t(pp_0)
    
    #robust fitting
    # TODO: What is this 2 and do I need to worry about it?
    res2 <- irls(fc_0,
                 w0_0,
                 probes,
                 ag = 2,
                 tol = 1e-3,
                 maxiter = 50)
    irlsLogStrings = res2[[4]]
    
    fp <- res2[[1]] #these are probe fitnesses fp
    #since fp is determined up to an additive constant, set the constant by
    #requiring that mean(fp[1:3]) = 0 (the null probes have zero fitness)
    # ab: so, the three null probes are first in the fp array because their names are, or start with, zero, so
    # they get put at the beginning by an alphabetical sort!
    # TODO: got to refactor out this 3
    mnull <- mean(fp[1:3])
    
    fp <- fp - mnull
    fc <-
      fc - mnull * 2 # TODO: What is this 2 and do I need to worry about it?
    # ah, I think I see where that 2 comes from: fp is at the probe level, but fc is at the *construct* level, and
    # each *construct* has *two* probes in it; we're subtracting out the value we'd expect for a construct made out of
    # *two* null probes.  So, no, I don't have to refactor out this two, because it is about the construct format
    # ("dual" crispr) not the number of replicates.
    
    #now shift a's according to the relation ac<-mx-mt*fc # ab: um, this comment isn't adjacent to any code! And wasn't the freedom removed from the ac's back in the fit_ac_fc function?
    #end removing freedom
    #a, fc, and fp are fully set
    
    rank_p <- rep(0, nprobes)
    names(rank_p) <- probes
    #find best probes
    fp12 <- fp
    i <- 1 #do null construct
    # TODO: refactor 3s
    rank_p[(i - 1) * 3 + 1:3] <-
      rank(abs(fp12[(i - 1) * 3 + 1:3])) #looking for the worst
    # eee--for the above expression, i is fixed at one.  Thus, the above expression simplifies to:
    # rank_p[(1-1)*3+1:3]<-rank(abs(fp12[(1-1)*3+1:3])) or
    # rank_p[1:3]<-rank(abs(fp12[1:3]))
    
    # TODO: refactor 3s (for # probes per gene)
    # The 2 here *isn't* actually number of replicates; it is just that we did first "gene"--really the null--
    # separately above.  The one difference between the equation used for the null and the one use for the
    # real genes is that for the null, abs is positive, whereas in the equation for the real genes, we are
    # examining the *negative* of the abs value.  This is because, for the null gene, we expect the real value
    # of fp for each of its probes to be zero--no effect on fitness.  Rank assigns ranks values in ascending order, so we
    # make them all positive using abs and say the one closest to zero is the best, the second closest to zero is
    # the second best, and so forth.  For example, if abs(fp12[1:3]) is 0.290, 0.009, 1.48, then
    # rank(abs(fp12[1:3])) = 2, 1, 3 since the first value in abs(fp12[1:3]) has a rank of 2 (second-closest to zero),
    # the econd value in abs(fp12[1:3]) has a rank of 1 (closest to zero), and so forth.
    # However, for the real genes, we ... assume? ... that the probe with the fp
    # value *farthest* from zero is the best probe.  So we make all the fp values positive using abs, then make
    # those abs values all negative.  Since rank assigns rank values in ascending order, the one with the most negative value
    # (i.e., the largest absolute value of fp) will be ranked first, etc.
    for (i in 2:n) {
      # for each gene, except the first one (which is really the null)
      rank_p[(i - 1) * 3 + 1:3] <-
        rank(-abs(fp12[(i - 1) * 3 + 1:3])) #looking for the best
    }
    
    #TODO: refactor 3, for number of probes per gene
    p_rank <-
      3 - rank_p # p_rank has the probes for each gene in reverse order from the way they are in rank_p; thus, the *best* probe in p_rank has the value 2 (3-1) whereas the worst has the value 0 (3-3), while in rank_p, the *best* probe has value 1 and the worst has value 3; it is used only to make the probe rank output file.
    
    wpi1 <- matrix(0, nrow = nprobes, ncol = nprobes)
    # for each pair of probes
    for (i in 1:nprobes) {
      for (j in 1:nprobes) {
        #TODO: refactor 3s
        #err ... subtracting 3 from rank_p will always either give a negative or zero.  But i guess it doesn't matter, as we're subtracting 3 from the rank of *both* probes and then multiplying, so the result will always be either positive or zero
        wpi1[i, j] <-
          (rank_p[i] - 3) * (rank_p[j] - 3) # product of reversed rank (i.e., best probe has biggest number) for the two probes in this probe pair
      }
    }
    
    f <- rep(0, n)
    names(f) <- genes
    for (i in 1:n) {
      # for each gene
      # TODO: refactor 3s
      w1 <-
        (rank_p[(i - 1) * 3 + 1:3] - 3) ^ 2 #ansatz for weights # so, subtract 3 from the rank of each of the probes from this gene, then square (squaring since value will be neg or zero--see above?)
      f[i] <-
        sum(w1 * fp[(i - 1) * 3 + 1:3]) / sum(w1) #weighted mean # multiply fp for each probe for this gene by the anzatz weight for that probe, then sum across all the probes for this gene, and divide that sum by the sum of the anzatz weights for all probes for this gene
      
    }
    
    pi1 <- res2[[3]] #raw pi-scores per construct
    
    mean_pi1 <- matrix(0, nrow = n, ncol = n)
    # oh, ffs, again with the for (i in 1:(n-1)) { ... for (j in (i+1):n) { construct
    for (i in 1:(n - 1)) {
      # for each first gene
      # TODO: refactor 3s
      ixi <- 3 * (i - 1) + 1:3 # get range of probe indices for first gene
      for (j in (i + 1):n) {
        # for each second gene not already covered
        ixj <- 3 * (j - 1) + 1:3 # get range of probe indices for second gene
        # reminder: w0_0 is the "physical goodness" of each construct.  Seems to be integer booleans (i.e., 1 for true, 0 for false)
        expressed1 <-
          w0_0[ixi, ixj] > 0 #define expressed probe pairs # ab: for this pair of genes
        
        # TODO: Not really sure what these two lines are doing; need to step through, understand
        # sum(wpi1[ixi,ixj][expressed1]), understand what kind of matrix calculations are happening
        local_w1 <-
          wpi1[ixi, ixj] / sum(wpi1[ixi, ixj][expressed1]) * sum(expressed1)
        # wait, local_w1 doesn't seem to actually be *used* anywhere ...?
        
        # the denominator here is clearly just making sure we never get a divide-by-zero error ...
        # looks like we're getting mean of the weighted pi values across all expressed probes for this gene
        mean_pi1[i, j] <-
          sum((pi1[ixi, ixj] * wpi1[ixi, ixj])[expressed1]) / max(small, sum(wpi1[ixi, ixj][expressed1]))
      }
    }
    
    uutri <- upper.tri(mean_pi1)
    uutri[1, ] <- FALSE #remove top line, 0 # why?
    zi1 <- mean_pi1[uutri]
    #zi <- zi1
    npi <- length(zi1)
    
    mmm <- length(fp) # number of probes
    pi_iter <- matrix(0, nrow = npi, ncol = niter)
    fp_iter <- matrix(0, nrow = mmm, ncol = niter)
    f_iter <- matrix(0, nrow = n, ncol = niter)
    
    utri <- upper.tri(fc_0)
    ntri <- sum(utri)
    #ppi_iter <- matrix(0, nrow = ntri, ncol = niter)
    
    for (iter in 1:niter) {
      irlsLogStrings = append(irlsLogStrings, capture.output(cat("\n", iter, "\n")))
      
      fc_1 <-
        matrix(0, nrow = nprobes, ncol = nprobes) # same starting value as fc_0<-matrix(0,nrow=nprobes,ncol=nprobes)
      #TODO: need to seed this to make it deterministic for testing
      set.seed(iter) # TODO: remove after testing.  This is so each time through loop has different--but still deteriministic--random numbers
      fc0 <-
        fc_0[utri] + rnorm(ntri, sd = sdfc_0[utri]) # ok, previous fc0 was fc0<-fc ; now we're adding a random normal to each fc, where each normal variable's mean is ?the sum of all fcs in the upper triangle?, and its stddev is the stddev of the fc of the analogous construct?
      pp0 <-
        pp_0[utri] # gets the posterior probability of the fcs for constructs in the upper triangle
      set.seed(iter) # TODO: remove after testing.  This is so each time through loop has different--but still deteriministic--random numbers
      draw <-
        ifelse(runif(ntri) < pp0, 1, 0) #TODO: need to seed this to make it deterministic for testing # draw is 1 if the random posterior probability is less than the calculated posterior probability, zero otherwise
      fc_1[utri] <-
        fc0 * draw # ok, multiplying by draw will set to zero every value in fc_1 where the posterior probability of the real fc was not more than you'd expect by chance
      fc_1 <- fc_1 + t(fc_1) # make fc_1 symmetric
      
      
      # TODO: repeat of code to run irls, retrieve outputs, constrain fp0
      #robust fitting
      # TODO: Do I need to worry about this 2?
      res2 <- irls(fc_1,
                   w0_0,
                   probes,
                   ag = 2,
                   tol = 1e-3,
                   maxiter = 50)
      irlsLogStrings = append(irlsLogStrings, res2[[4]])
      
      fp0 <- res2[[1]] #these are probe fitnesses fp
      
      #since fp is determined up to an additive constant, set the constant by
      #requiring that mean(fp[1:3]) = 0 (the null probes have zero fitness
      # TODO: refactor 3
      mnull <- mean(fp0[1:3])
      fp0 <- fp0 - mnull
      # end repeat of code
      
      # TODO: near-repeat of code to set weighted mean f above; now have one f per iteration
      for (i in 1:n) {
        # TODO: refactor 3s
        w1 <- (rank_p[(i - 1) * 3 + 1:3] - 3) ^ 2 #ansatz for weights
        f_iter[i, iter] <-
          sum(w1 * fp0[(i - 1) * 3 + 1:3]) / sum(w1) #weighted mean
      }
      # end near-repeat of code
      
      pi1 <- res2[[3]] #raw pi-scores per construct # repeat code
      pi_scrambled <-
        pi1 # why put pi1 into pi_scrambled, then use it without changing it in any way? It *isn't* scrambled, it is just pi1
      
      # TODO: repeat of above code to set mean_pi1
      mean_pi1 <- matrix(0, nrow = n, ncol = n)
      for (i in 1:(n - 1)) {
        # TODO: refactor 3s
        ixi <- 3 * (i - 1) + 1:3
        for (j in (i + 1):n) {
          ixj <- 3 * (j - 1) + 1:3
          expressed1 <- w0_0[ixi, ixj] > 0 #define expressed probe pairs
          local_w1 <-
            wpi1[ixi, ixj] / sum(wpi1[ixi, ixj][expressed1]) * sum(expressed1)
          
          mean_pi1[i, j] <-
            sum((pi_scrambled[ixi, ixj] * wpi1[ixi, ixj])[expressed1]) / max(small, sum(wpi1[ixi, ixj][expressed1]))
        }
      }
      zi1 <- mean_pi1[uutri]
      # end repeat of code
      
      pi_iter[, iter] <- zi1
      fp_iter[, iter] <- fp0
    }
    
    #f_mean<-apply(f_iter,1,mean) # one fmean for each gene # wait, this value is never used!  What is output in the single-gene fitness file is f, not f_mean, although the stddev output there is f_sd (see right below)
    f_sd <- apply(f_iter, 1, sd) # output in single gene fitness file
    
    #fp_mean<-apply(fp_iter,1,mean) # never used
    #fp_sd<-apply(fp_iter,1,sd) # never used
    
    pi_mean <-
      apply(pi_iter, 1, mean) # used in several plots, output in pi score file
    pi_sd <- apply(pi_iter, 1, sd) # output in pi score file
    pi_iter_null <-
      pi_iter - pi_mean #used directly below, and in prepping for pi score output file
    pi_null <-
      c(pi_iter_null, -pi_iter_null) # used directly below, and in histogram of pi scores
    return(
      list(
        pi_mean = pi_mean,
        pi_null = pi_null,
        pi_sd = pi_sd,
        pi_iter_null = pi_iter_null,
        uutri = uutri,
        fc = fc,
        f = f,
        p_rank = p_rank,
        f_sd = f_sd,
        irlsLogStrings = irlsLogStrings
      )
    )
  }