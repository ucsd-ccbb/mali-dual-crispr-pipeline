# input params
niter = 2 #TODO: put back to 1000 once determine that testing framework is sound
g_cran_repo_url = 'http://cran.us.r-project.org'
time = c(3,14,20,28)
project = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr"
# TODO: Abundance-level hardcode
ab0 <- c(-19. , -18.5, -18.5, -19. , -19. , -19. , -19. , -19.)
# end input params

library(qvalue)

X<-read.table(input_filename,sep="\t",header=TRUE)

small<-1e-6
nt<-length(time) #nt >= 2 # number of timepoints
mt<-mean(time)
vt<-Var(time)


tempPrepDataOutput = prepData(input_filename, nt)
nn=tempPrepDataOutput$nn
probes=tempPrepDataOutput$probes
nprobes=tempPrepDataOutput$nprobes
genes=tempPrepDataOutput$genes
n=tempPrepDataOutput$n
pA_pB=tempPrepDataOutput$pA_pB
y=tempPrepDataOutput$y # y is df of log2 frequencies (rows = constructs, cols = timept+replicate, no info columns)

tempXsAbsAndBads = getXsAbsAndBads(y, nt, ab0)
x1 = tempXsAbsAndBads$x1 # x1 is log2 frequencies for the 1st replicate of all timepts
ab1 = tempXsAbsAndBads$ab1 # ab1 is abundance thresholds for all 1st replicates
bad1 = tempXsAbsAndBads$bad1 
x2 = tempXsAbsAndBads$x2 # x2 is log2 frequencies for the 2nd replicate of all timepts
ab2 = tempXsAbsAndBads$ab2 # ab2 is abundance thresholds for all 2nd replicates
bad2 = tempXsAbsAndBads$bad2
rownames(x1)<-pA_pB 
rownames(x2)<-pA_pB

# changed from print to output to (temporary) global variable to be used in testing
# TODO: when refactoring code, change these values to output of function, then update test
gOutputBad1Bad2 = c(sum(bad1), sum(bad2))
# changed from print to global variable for testing
# TODO: refactor to be output of function
gSumBad1Bad2 = sum(bad1 & bad2)

# TODO: "2 replicate" id is baked right into this function definition.
resf<-fit_ac_fc(x1,ab1,x2,ab2) # I don't know what resf means ... maybe "resulting fit"?  Anyway, it is a list, which is unpacked below.

# TODO: 2-replicate assumption is baked into outputs of fit_ac_fc function
a1<-resf[[1]] # a1 is the initial condition (in log2 frequency) for each construct c for replicate 1
a2<-resf[[2]] # a2 is the initial condition (in log2 frequency) for each construct c for replicate 2
fc<-resf[[3]] # fc is the fitness of each construct c (calculated across both replicates)
sdfc<-resf[[4]] #standard error # Roman's comment says this is stderr, but I don't think it actually is--I think sdfc is the std deviation of the fitness of each construct c (calculated across both replicates)
p_t<-resf[[5]] #raw p-value from t-test # p_t is the raw p value of the fc of each construct c (calculated across both replicates)
lfdr_fc<-resf[[6]] #lfdr from p_t (Storey) # lfdr_fc is the local FDR of each construct (calculated across both replicates)
pp_fc<-1-lfdr_fc  # ab: I believe this is posterior probability of fc of each construct c (calculated across both replicates)
# Wikipedia: "the posterior probability of ... an uncertain proposition is the conditional probability that is assigned after the relevant evidence ... is taken into account".  Thus, higher equals more probable.
df<-resf[[7]] #degrees of freedom # df is the degrees of freedom of each construct c (calculated across both replicates)
allbad<-resf[[8]] #is TRUE when both experiments are bad (at most 1 good value) # allbad is a boolean value for each construct c that is true for all the constructs that lack at least 2 acceptable-abundance timepoints in BOTH experiments

# assign construct names to the newly-created ... arrays? ... of per construct info
names(fc)<-pA_pB
names(pp_fc)<-pA_pB
names(sdfc)<-pA_pB

# plot fc for all good constructs vs posterior probability of fc for all those same good constructs
plot(fc[!allbad],pp_fc[!allbad],pch=16,cex=0.2,xlab=expression(f["c"]),ylab="posterior probability")
plotRandomVsRealPosteriorProbs(nn, fc, pp_fc)

# TODO: I suspect that the "2^"s here are reversing a log2 transform rather than
# hardcoding the # of replicates, so they are probably ok.  However, the repeat of the
# hist 2 times is clearly for the 2 replicates and needs to be refactored.
hist(2^a1,breaks=1000,xlab="relative abundance",main="replicate 1, time 0")
hist(2^a2,breaks=1000,xlab="relative abundance",main="replicate 2, time 0")

tempRobustFittingOutput = doRobustFitting(nn, pA_pB, allbad, fc, nprobes, probes, n, pp_fc, sdfc)
pi_mean = tempRobustFittingOutput$pi_mean
pi_null = tempRobustFittingOutput$pi_null
pi_sd = tempRobustFittingOutput$pi_sd
pi_iter_null = tempRobustFittingOutput$pi_iter_null
uutri = tempRobustFittingOutput$uutri
fc = tempRobustFittingOutput$fc
f = tempRobustFittingOutput$f
p_rank = tempRobustFittingOutput$p_rank
f_sd = tempRobustFittingOutput$f_sd
gIrlsLogStrings = tempRobustFittingOutput$irlsLogStrings

rge<-range(pi_mean)
# TODO: where does 0.002 in break seq come from?  Should we set in params?
h<-hist(pi_mean,breaks=seq(rge[1]-0.002,rge[2]+0.002,by=0.002),main="",xlab=expression(pi["gg'"]),col="grey80",border=FALSE,probability=TRUE)
d<-density(pi_null,bw=0.002)
lines(d,col="black")
rug(pi_mean)
box()

gFdrsList = calcFdrsLeftAndRight(pi_mean, pi_null)
plotFdrsLeftAndRight(pi_mean, pi_null, gFdrsList$fdr_left, gFdrsList$fdr_right)
writePosteriorProbabilityFile(n, genes, pi_mean, pi_sd, pi_iter_null, uutri, f, gFdrsList$fdr_left, gFdrsList$fdr_right)
plotFitToPdf(x1,a1,fc,ab1,x2,a2,ab2)
writeConstructFitnessesFile(pA_pB, fc, sdfc)
plotConstructFitnessesHist(fc, allbad)
plot_scatterplots()
writeProbeRanksFile(probes,p_rank)
writeSingleGeneFitnessesFile(genes, f, f_sd)
plotSingleGeneFitnesses(f, genes)