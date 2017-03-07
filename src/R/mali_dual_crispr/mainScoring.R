scoreDualCrisprScreen <- function(input_filename, niter, time, project, ab0) {
  numTimepoints <- length(time) #numTimepoints >= 2 # number of timepoints
  
  # TODO: Determine if still need to do data prep in R or if all needed is covered by what I put in scoring prep output file :)
  tempPrepDataOutput = prepData(input_filename, numTimepoints)
  numDiffGeneConstructs = tempPrepDataOutput$nn
  uniqueAlphabeticalProbesList = tempPrepDataOutput$probes
  numUniqueProbes = tempPrepDataOutput$nprobes
  uniqueAlphabeticalGenesList = tempPrepDataOutput$genes
  numUniqueGenes = tempPrepDataOutput$n
  constructNamesList = tempPrepDataOutput$pA_pB
  constructBySampleLog2FreqsMatrix = tempPrepDataOutput$y # constructBySampleLog2FreqsMatrix is df of log2 frequencies (rows = constructs, cols = timept+replicate, no info columns)
  
  tempXsAbsAndBads = getXsAbsAndBads(constructBySampleLog2FreqsMatrix, numTimepoints, ab0)
  rep1constructBySampleLog2FreqsMatrix = tempXsAbsAndBads$x1 # rep1constructBySampleLog2FreqsMatrix is log2 frequencies for the 1st replicate of all timepts
  rep1abundanceThresholds = tempXsAbsAndBads$ab1 # rep1abundanceThresholds is abundance thresholds for all 1st replicates
  rep1LtTwoTimeptsAboveAbundanceThresholdMask = tempXsAbsAndBads$bad1
  rep2constructBySampleLog2FreqsMatrix = tempXsAbsAndBads$x2 # rep2constructBySampleLog2FreqsMatrix is log2 frequencies for the 2nd replicate of all timepts
  rep2abundanceThresholds = tempXsAbsAndBads$ab2 # rep2abundanceThresholds is abundance thresholds for all 2nd replicates
  rep2LtTwoTimeptsAboveAbundanceThresholdMask = tempXsAbsAndBads$bad2
  rownames(rep1constructBySampleLog2FreqsMatrix) <- constructNamesList
  rownames(rep2constructBySampleLog2FreqsMatrix) <- constructNamesList
  
  # TODO: "2 replicate" id is baked right into this function definition.
  resf <-
    fit_ac_fc(rep1constructBySampleLog2FreqsMatrix, rep1abundanceThresholds, rep2constructBySampleLog2FreqsMatrix, rep2abundanceThresholds, numTimepoints) # I don't know what resf means ... maybe "resulting fit"?  Anyway, it is a list, which is unpacked below.
  
  # TODO: 2-replicate assumption is baked into outputs of fit_ac_fc function
  a1 <-
    resf[[1]] # a1 is the initial condition (in log2 frequency) for each construct c for replicate 1
  a2 <-
    resf[[2]] # a2 is the initial condition (in log2 frequency) for each construct c for replicate 2
  fc <-
    resf[[3]] # fc is the fitness of each construct c (calculated across both replicates)
  sdfc <-
    resf[[4]] #standard error # Roman's comment says this is stderr, but I don't think it actually is--I think sdfc is the std deviation of the fitness of each construct c (calculated across both replicates)
  lfdr_fc <-
    resf[[6]] #lfdr from p_t (Storey) # lfdr_fc is the local FDR of each construct (calculated across both replicates)
  pp_fc <-
    1 - lfdr_fc  # ab: I believe this is posterior probability of fc of each construct c (calculated across both replicates)
  # Wikipedia: "the posterior probability of ... an uncertain proposition is the conditional probability that is assigned after the relevant evidence ... is taken into account".  Thus, higher equals more probable.
  allbad <-
    resf[[8]] #is TRUE when both experiments are bad (at most 1 good value) # allbad is a boolean value for each construct c that is true for all the constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments
  
  # assign construct names to the newly-created ... arrays? ... of per construct info
  names(fc) <- constructNamesList
  names(pp_fc) <- constructNamesList
  names(sdfc) <- constructNamesList
  
  # plot fc for all good constructs vs posterior probability of fc for all those same good constructs
  plot(
    fc[!allbad],
    pp_fc[!allbad],
    pch = 16,
    cex = 0.2,
    xlab = expression(f["c"]),
    ylab = "posterior probability"
  )
  plotRandomVsRealPosteriorProbs(numDiffGeneConstructs, fc, pp_fc, allbad)
  
  # TODO: I suspect that the "2^"s here are reversing a log2 transform rather than
  # hardcoding the # of replicates, so they are probably ok.  However, the repeat of the
  # hist 2 times is clearly for the 2 replicates and needs to be refactored.
  hist(2 ^ a1,
       breaks = 1000,
       xlab = "relative abundance",
       main = "replicate 1, time 0")
  hist(2 ^ a2,
       breaks = 1000,
       xlab = "relative abundance",
       main = "replicate 2, time 0")
  
  tempRobustFittingOutput = doRobustFitting(numDiffGeneConstructs, constructNamesList, allbad, fc, numUniqueProbes, uniqueAlphabeticalProbesList, numUniqueGenes, pp_fc, sdfc, uniqueAlphabeticalGenesList, niter)
  pi_mean = tempRobustFittingOutput$pi_mean
  pi_null = tempRobustFittingOutput$pi_null
  pi_sd = tempRobustFittingOutput$pi_sd
  pi_iter_null = tempRobustFittingOutput$pi_iter_null
  uutri = tempRobustFittingOutput$uutri
  fc = tempRobustFittingOutput$fc
  f = tempRobustFittingOutput$f
  p_rank = tempRobustFittingOutput$p_rank
  f_sd = tempRobustFittingOutput$f_sd
  
  rge <- range(pi_mean)
  # TODO: where does 0.002 in break seq come from?  Should we set in params?
  hist(
      pi_mean,
      breaks = seq(rge[1] - 0.002, rge[2] + 0.002, by = 0.002),
      main = "",
      xlab = expression(pi["gg'"]),
      col = "grey80",
      border = FALSE,
      probability = TRUE
    )
  d <- density(pi_null, bw = 0.002)
  lines(d, col = "black")
  rug(pi_mean)
  box()
  
  fdrsList = calcFdrsLeftAndRight(pi_mean, pi_null)
  plotFdrsLeftAndRight(pi_mean, pi_null, fdrsList$fdr_left, fdrsList$fdr_right)
  writePosteriorProbabilityFile(project,
    numUniqueGenes,
    uniqueAlphabeticalGenesList,
    pi_mean,
    pi_sd,
    pi_iter_null,
    uutri,
    f,
    fdrsList$fdr_left,
    fdrsList$fdr_right
  )
  plotFitToPdf(project, rep1constructBySampleLog2FreqsMatrix, a1, fc, rep1abundanceThresholds, rep2constructBySampleLog2FreqsMatrix, a2, rep2abundanceThresholds, numTimepoints, constructNamesList)
  writeConstructFitnessesFile(project, constructNamesList, fc, sdfc)
  plotConstructFitnessesHist(fc, allbad)
  plot_scatterplots(numDiffGeneConstructs, numTimepoints, a1, a2, fc, time, rep1LtTwoTimeptsAboveAbundanceThresholdMask, rep2LtTwoTimeptsAboveAbundanceThresholdMask, rep1constructBySampleLog2FreqsMatrix, rep2constructBySampleLog2FreqsMatrix)
  writeProbeRanksFile(project, uniqueAlphabeticalProbesList, p_rank)
  writeSingleGeneFitnessesFile(project, uniqueAlphabeticalGenesList, f, f_sd)
  plotSingleGeneFitnesses(f, uniqueAlphabeticalGenesList)
  
  return(
    list(
      rep1LtTwoTimeptsAboveAbundanceThresholdMask = rep1LtTwoTimeptsAboveAbundanceThresholdMask,
      rep2LtTwoTimeptsAboveAbundanceThresholdMask = rep2LtTwoTimeptsAboveAbundanceThresholdMask,
      irlsLogStrings = tempRobustFittingOutput$irlsLogStrings
    )
  )
}