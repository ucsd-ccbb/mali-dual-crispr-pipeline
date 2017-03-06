# input params
niter = 2 #TODO: put back to 1000 once determine that testing framework is sound
g_cran_repo_url = 'http://cran.us.r-project.org'
input_filename = "/Temp/A549_CV4_counts_w_everything.txt"
time = c(3,14,20,28)
project = "A549_CV4_3-14-21-28_NA_combined_simple-null-w-lfdr"
# end input params

library(qvalue)

# sourcing shared functions
source('sharedRfunctions.R')
# end sourcing

X<-read.table(input_filename,sep="\t",header=TRUE)

small<-1e-6
nt<-length(time) #nt >= 2 # number of timepoints
mt<-mean(time)
vt<-Var(time)

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

# TODO: Abundance-level hardcode
ab0 <- c(-19. , -18.5, -18.5, -19. , -19. , -19. , -19. , -19.)

# So ... we're building 3 matrices of n genes by n genes.
# The first one just 
# Again, here 1 & 2 hardcoding acceptable as means 1st and 2nd of DUAL crispr construct
g1names<-matrix("",ncol=n,nrow=n) # n = number of genes
g2names<-matrix("",ncol=n,nrow=n)
ggnames<-matrix("",ncol=n,nrow=n)
for (i in 1:(n-1)) {
  for (j in (i+1):n) {
    g1names[i,j]<-genes[i]
    g2names[i,j]<-genes[j]
    # TODO: Separator should be set in params rather than hardcoded
    ggnames[i,j]<-paste(genes[i],"_",genes[j],sep="")
    ggnames[j,i]<-ggnames[i,j]
  }
}
# OK, everything above this comment actually has zero to do with finding the # of constructs below the
# abundance thresholds--it is just setting up some data structures with gene labels in them for use (much) later
# in the pi score file output.


# y is df of log2 frequencies (rows = constructs, cols = timept+replicate, no info columns)
# ab0 is the vector of log2 frequency abundance thresholds, e.g.,
# array([-19. , -18.5, -18.5, -19. , -19. , -19. , -19. , -19. ])

# for replicate 1
# "2"s here and line below look like replicate hardcoding; refactor

# Ok, so this is getting the log2 frequencies for the 1st replicate of all timepts
x1<-y[,seq(1,2*nt,by=2)]
# This is getting the abundance thresholds for all of these 1st replicates
ab1<-ab0[seq(1,2*nt,by=2)]
# TODO: 1st "2" here is ok (designates row vs col) ... but what is source of second?
bad1<-apply(t(x1)>ab1,2,sum)<2
print(sum(bad1))

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
print(sum(bad2))

# the printed numbers are *per-replicate* ... I guess we wanted them that way and that is why we did 
# the above process separately for each replicate?

# x1 is log2 frequencies for the 1st replicate of all timepts
# x2 is log2 frequencies for the 2nd replicate of all timepts
# ab1 is abundance thresholds for all 1st replicates
# ab2 is abundance thresholds for all 2nd replicates
# after this cell, these 4 are used again one more time considerably later--in the call to plot_fit that makes
# Log2 Frequency vs Time Plots for Constructs with Large Fitness File Output

rownames(x1)<-pA_pB 
rownames(x2)<-pA_pB

# TODO: ok, I don't know what this function does, but since ab1 and ab2 are calculated for replicates 1 and 2,
# I think that probably the "2 replicate" id is baked right into this function definition.
resf<-fit_ac_fc(x1,ab1,x2,ab2)
# I don't know what resf means ... maybe "resulting fit"?  Anyway, it is a list, which is unpacked below.

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

# plot fc for all good constructs vs posterior probability of fc for all those same good constructs
plot(fc[!allbad],pp_fc[!allbad],pch=16,cex=0.2,xlab=expression(f["c"]),ylab="posterior probability")

# nn = number of *good* constructs
# r = nn random numbers pulled from a uniform distribution with (default) min of 0 and max of 1 
set.seed(1) # TODO: remove after testing.  
r<-runif(nn) #TODO: need to seed this to make it deterministic for testing
fr<-fc[r<pp_fc] # errr ... fc for all the constructs where the random ... posterior probability? ... is less than the actual calculated posterior probability for this construct?
# thus, fr is the distribution of fcs for constructs whose true pp of happening is greater than their random pp of happening--roughly, fcs of constructs whose fc values are more likely than chance
rge<-range(fc) # range of fcs across all constructs

# TODO: what determines the 0.001 for the break sequence?  Need to set in params?
plotOverlappingHist(fc[!allbad],fr,breaks=seq(rge[1]-.001,rge[2]+0.001,by=0.001),xlab=expression(f["c"]),ylab="Frequency")

# assign construct names to the newly-created ... arrays? ... of per construct info
names(fc)<-pA_pB
names(pp_fc)<-pA_pB
names(sdfc)<-pA_pB
# TODO: I suspect that the "2^"s here are reversing a log2 transform rather than
# hardcoding the # of replicates, so they are probably ok.  However, the repeat of the
# hist 2 times is clearly for the 2 replicates and needs to be refactored.
hist(2^a1,breaks=1000,xlab="relative abundance",main="replicate 1, time 0")
hist(2^a2,breaks=1000,xlab="relative abundance",main="replicate 2, time 0")

# nn = number of *good* constructs
# whatever u is, it is 0 for all bad and/or non-existent constructs and, at least at the beginning, one for all good constructs
u1<-rep(0,nn)
names(u1)<-pA_pB
u1[!allbad]<-1  #all other weights set to 1

# argh, heaven preserve me from variables fc0 and fc_0 that mean different things :(
fc0<-fc 

# just sets expected matrix size for all these--num probes by num probes, w/defaults of zero, and row/col names
fc_0<-matrix(0,nrow=nprobes,ncol=nprobes)
sdfc_0<-matrix(0,nrow=nprobes,ncol=nprobes)
w0_0<-matrix(0,nrow=nprobes,ncol=nprobes)
pp_0<-matrix(0,nrow=nprobes,ncol=nprobes)
rownames(fc_0)<-probes
colnames(fc_0)<-probes
rownames(sdfc_0)<-probes
colnames(sdfc_0)<-probes
rownames(w0_0)<-probes
colnames(w0_0)<-probes
rownames(pp_0)<-probes
colnames(pp_0)<-probes


# for 1 to number of genes - 1
for (i in 1:(n-1)) {
  # TODO: Ok, where did these "3"s come from?  From the expected number of probes per gene?
  # If so, that's *got* to be refactored out ...
  for (k in 1:3) {
    # index of current probe for current gene in probes array
    iprobe<-(i-1)*3+k
    # probes is entire probe name set in alphabetical order
    iprobe_name<-probes[iprobe]
    
    # note: for (i in 1:(n-1)) {
    #for (j in (i+1):n) {
    # nested loops here are the same as used above to generate ggnames...
    # for all possible second genes that come after the current first gene (i) 
    for (j in (i+1):n) {
      # for all the constructs of the second gene
      for (l in 1:3) {
        jprobe<-(j-1)*3+l
        jprobe_name<-probes[jprobe] # as for first gene, get index and then name of second gene
        
        #generate the construct name
        construct<-paste(iprobe_name,"_",jprobe_name,sep="")
        w0_0[iprobe_name,jprobe_name]<-u1[construct] #initial weights. non-existent pairs will have w0=0.  #ab: Roman's comment here says this is initial weights, but his comment in the irls function says that this is actually the "physical goodness" of each construct.  Seems to be integer booleans (i.e., 1 for true, 0 for false) 
        pp_0[iprobe_name,jprobe_name]<-pp_fc[construct] #initial weights. non-existent pairs will have w0=0 # ab: I don't think this comment is right either.  pp_fc are the posterior probabilities for the fc of each construct
        fc_0[iprobe_name,jprobe_name]<-fc0[construct]
        sdfc_0[iprobe_name,jprobe_name]<-sdfc[construct]
      }
    }
  }
}
w0_0<-w0_0+t(w0_0) #make symmetric
fc_0<-fc_0+t(fc_0)
pp_0<-pp_0+t(pp_0)

#robust fitting
# TODO: What is this 2 and do I need to worry about it?
res2<-irls(fc_0,w0_0,probes,ag=2,tol=1e-3,maxit=50)

fp<-res2[[1]] #these are probe fitnesses fp
#since fp is determined up to an additive constant, set the constant by 
#requiring that mean(fp[1:3]) = 0 (the null probes have zero fitness)
# ab: so, the three null probes are first in the fp array because their names are, or start with, zero, so
# they get put at the beginning by an alphabetical sort!
# TODO: got to refactor out this 3
mnull<-mean(fp[1:3])

fp<-fp-mnull
fc<-fc-mnull*2 # TODO: What is this 2 and do I need to worry about it?
# ah, I think I see where that 2 comes from: fp is at the probe level, but fc is at the *construct* level, and
# each *construct* has *two* probes in it; we're subtracting out the value we'd expect for a construct made out of
# *two* null probes.  So, no, I don't have to refactor out this two, because it is about the construct format 
# ("dual" crispr) not the number of replicates.

#now shift a's according to the relation ac<-mx-mt*fc # ab: um, this comment isn't adjacent to any code! And wasn't the freedom removed from the ac's back in the fit_ac_fc function?
#end removing freedom
#a, fc, and fp are fully set

rank_p<-rep(0,nprobes)
names(rank_p)<-probes
#find best probes
fp12<-fp
i<-1 #do null construct
# TODO: refactor 3s
rank_p[(i-1)*3+1:3]<-rank(abs(fp12[(i-1)*3+1:3])) #looking for the worst
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
for (i in 2:n) { # for each gene, except the first one (which is really the null)
  rank_p[(i-1)*3+1:3]<-rank(-abs(fp12[(i-1)*3+1:3])) #looking for the best
}

#TODO: refactor 3, for number of probes per gene
p_rank<-3-rank_p # p_rank has the probes for each gene in reverse order from the way they are in rank_p; thus, the *best* probe in p_rank has the value 2 (3-1) whereas the worst has the value 0 (3-3), while in rank_p, the *best* probe has value 1 and the worst has value 3; it is used only to make the probe rank output file.

wpi1<-matrix(0,nrow=nprobes,ncol=nprobes)
# for each pair of probes
for (i in 1:nprobes) {
  for (j in 1:nprobes) {
    #TODO: refactor 3s
    #err ... subtracting 3 from rank_p will always either give a negative or zero.  But i guess it doesn't matter, as we're subtracting 3 from the rank of *both* probes and then multiplying, so the result will always be either positive or zero
    wpi1[i,j]<-(rank_p[i]-3)*(rank_p[j]-3) # product of reversed rank (i.e., best probe has biggest number) for the two probes in this probe pair
  }
}

f<-rep(0,n)
names(f)<-genes
for (i in 1:n) { # for each gene
  # TODO: refactor 3s
  w1<-(rank_p[(i-1)*3+1:3]-3)^2 #ansatz for weights # so, subtract 3 from the rank of each of the probes from this gene, then square (squaring since value will be neg or zero--see above?)
  f[i]<-sum(w1*fp[(i-1)*3+1:3])/sum(w1) #weighted mean # multiply fp for each probe for this gene by the anzatz weight for that probe, then sum across all the probes for this gene, and divide that sum by the sum of the anzatz weights for all probes for this gene
  
}
fmean<-f # fmean is array of zero values, one for each gene

pi1<-res2[[3]] #raw pi-scores per construct

mean_pi1<-matrix(0,nrow=n,ncol=n)
# oh, ffs, again with the for (i in 1:(n-1)) { ... for (j in (i+1):n) { construct
for (i in 1:(n-1)) { # for each first gene 
  # TODO: refactor 3s
  ixi<-3*(i-1)+1:3 # get range of probe indices for first gene 
  for (j in (i+1):n) { # for each second gene not already covered
    ixj<-3*(j-1)+1:3 # get range of probe indices for second gene 
    # reminder: w0_0 is the "physical goodness" of each construct.  Seems to be integer booleans (i.e., 1 for true, 0 for false)
    expressed1<-w0_0[ixi,ixj]>0 #define expressed probe pairs # ab: for this pair of genes
    
    # TODO: Not really sure what these two lines are doing; need to step through, understand 
    # sum(wpi1[ixi,ixj][expressed1]), understand what kind of matrix calculations are happening
    local_w1<-wpi1[ixi,ixj]/sum(wpi1[ixi,ixj][expressed1])*sum(expressed1)
    # wait, local_w1 doesn't seem to actually be *used* anywhere ...?
    
    # the denominator here is clearly just making sure we never get a divide-by-zero error ...
    # looks like we're getting mean of the weighted pi values across all expressed probes for this gene
    mean_pi1[i,j]<-sum((pi1[ixi,ixj]*wpi1[ixi,ixj])[expressed1])/max(small,sum(wpi1[ixi,ixj][expressed1]))
  }
}

uutri<-upper.tri(mean_pi1)
uutri[1,]<-FALSE #remove top line, 0 # why?
zi1<-mean_pi1[uutri]
zi<-zi1
npi<-length(zi1)

mmm<-length(fp) # number of probes
pi_iter<-matrix(0,nrow=npi,ncol=niter)
fp_iter<-matrix(0,nrow=mmm,ncol=niter)
f_iter<-matrix(0,nrow=n,ncol=niter)

utri<-upper.tri(fc_0)
ntri<-sum(utri)
ppi_iter<-matrix(0,nrow=ntri,ncol=niter)

for (iter in 1:niter) {
  cat("\n",iter,"\n")
  
  fc_1<-matrix(0,nrow=nprobes,ncol=nprobes) # same starting value as fc_0<-matrix(0,nrow=nprobes,ncol=nprobes)
  #TODO: need to seed this to make it deterministic for testing
  set.seed(iter) # TODO: remove after testing.  This is so each time through loop has different--but still deteriministic--random numbers
  fc0<-fc_0[utri]+rnorm(ntri,sd=sdfc_0[utri]) # ok, previous fc0 was fc0<-fc ; now we're adding a random normal to each fc, where each normal variable's mean is ?the sum of all fcs in the upper triangle?, and its stddev is the stddev of the fc of the analogous construct?
  pp0<-pp_0[utri] # gets the posterior probability of the fcs for constructs in the upper triangle
  set.seed(iter) # TODO: remove after testing.  This is so each time through loop has different--but still deteriministic--random numbers
  draw<-ifelse(runif(ntri)<pp0,1,0) #TODO: need to seed this to make it deterministic for testing # draw is 1 if the random posterior probability is less than the calculated posterior probability, zero otherwise
  fc_1[utri]<-fc0*draw # ok, multiplying by draw will set to zero every value in fc_1 where the posterior probability of the real fc was not more than you'd expect by chance
  fc_1<-fc_1+t(fc_1) # make fc_1 symmetric
  
  
  # TODO: repeat of code to run irls, retrieve outputs, constrain fp0
  #robust fitting
  # TODO: Do I need to worry about this 2?
  res2<-irls(fc_1,w0_0,probes,ag=2,tol=1e-3,maxit=50)
  
  fp0<-res2[[1]] #these are probe fitnesses fp
  
  #since fp is determined up to an additive constant, set the constant by 
  #requiring that mean(fp[1:3]) = 0 (the null probes have zero fitness
  # TODO: refactor 3
  mnull<-mean(fp0[1:3])
  fp0<-fp0-mnull
  # end repeat of code
  
  # TODO: near-repeat of code to set weighted mean f above; now have one f per iteration
  for (i in 1:n) {
    # TODO: refactor 3s
    w1<-(rank_p[(i-1)*3+1:3]-3)^2 #ansatz for weights
    f_iter[i,iter]<-sum(w1*fp0[(i-1)*3+1:3])/sum(w1) #weighted mean
  }
  # end near-repeat of code
  
  pi1<-res2[[3]] #raw pi-scores per construct # repeat code
  pi_scrambled<-pi1 # why put pi1 into pi_scrambled, then use it without changing it in any way? It *isn't* scrambled, it is just pi1
  
  # TODO: repeat of above code to set mean_pi1
  mean_pi1<-matrix(0,nrow=n,ncol=n)
  for (i in 1:(n-1)) {
    # TODO: refactor 3s
    ixi<-3*(i-1)+1:3
    for (j in (i+1):n) {
      ixj<-3*(j-1)+1:3
      expressed1<-w0_0[ixi,ixj]>0 #define expressed probe pairs
      local_w1<-wpi1[ixi,ixj]/sum(wpi1[ixi,ixj][expressed1])*sum(expressed1)
      
      mean_pi1[i,j]<-sum((pi_scrambled[ixi,ixj]*wpi1[ixi,ixj])[expressed1])/max(small,sum(wpi1[ixi,ixj][expressed1]))
    }
  }
  zi1<-mean_pi1[uutri]
  # end repeat of code
  
  pi_iter[,iter]<-zi1
  fp_iter[,iter]<-fp0
}

f_mean<-apply(f_iter,1,mean) # one fmean for each gene # wait, this value is never used!  What is output in the single-gene fitness file is f, not f_mean, although the stddev output there is f_sd (see right below)
f_sd<-apply(f_iter,1,sd) # output in single gene fitness file

fp_mean<-apply(fp_iter,1,mean) # never used
fp_sd<-apply(fp_iter,1,sd) # never used

pi_mean<-apply(pi_iter,1,mean) # used in several plots, output in pi score file
pi_sd<-apply(pi_iter,1,sd) # output in pi score file

pi_iter_null<-pi_iter-pi_mean #used directly below, then *redefined* (same way) and used (differently) in prepping for pi score output file

pi_null<-c(pi_iter_null,-pi_iter_null) # used directly below, and in histogram of pi scores
enull<-ecdf(pi_null) # used in plot of pi scores by fdr
emean<-ecdf(pi_mean) # used in plot of pi scores by fdr

rge<-range(pi_mean)
# TODO: where does 0.002 in break seq come from?  Should we set in params?
h<-hist(pi_mean,breaks=seq(rge[1]-0.002,rge[2]+0.002,by=0.002),main="",xlab=expression(pi["gg'"]),col="grey80",border=FALSE,probability=TRUE)
d<-density(pi_null,bw=0.002)
lines(d,col="black")
rug(pi_mean)
box()

fdr_left<-pmin(1,enull(pi_mean)/emean(pi_mean))
fdr_right<-pmin(1,(enull(-pi_mean))/(1-emean(pi_mean)))
plot(pi_mean,fdr_left,ylim=c(0,1),xlab=expression(pi["gg'"]),ylab="FDR, left") #left tail test
plot(pi_mean,fdr_right,ylim=c(0,1),xlab=expression(pi["gg'"]),ylab="FDR, right") #right tail test

z<-pi_mean/sd(pi_mean)

pi_iter_null<-pi_iter-pi_mean
abspi<-abs(pi_mean)
PP<-apply(abs(pi_iter_null)<abspi,1,mean)

oPP<-order(z)

#get names of these gene pairs
# These g1//g2, etc hardcodes are ok because represent DUAL crispr genes
names_of_g1<-g1names[uutri]
fg1<-f[names_of_g1]
names_of_g2<-g2names[uutri]
fg2<-f[names_of_g2]
fg12<-fg1+fg2
names_of_gg<-ggnames[uutri]

res<-data.frame(names_of_gg,names_of_g1,fg1,names_of_g2,fg2,fg12,pi_mean,pi_sd,PP,abspi,fdr_left,fdr_right,z)
# TODO: Let's refactor these colnames, file suffix, and separator somewhere easier to manage
colnames(res)<-c("gene_gene","geneA","fA","geneB","fB","fA+fB","pi","sd","PP","abs pi","FDR left","FDR right","z")
#TODO: remove format call here; added to aid in testing
write.table(format(res[oPP,], digits=14),file=paste(project,"_pi.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

pdf(paste(project,".pdf",sep=""))
par(pty="s",mfrow=c(1,1))
# TODO: I don't know what this is doing, but it looks like the method interface has 2 replicates baked right in ...
plot_fit(x1,a1,fc,ab1,x2,a2,fc,ab2,minfc=0.10)
dev.off()

resp<-data.frame(pA_pB,fc,sdfc)
# TODO: refactor file suffix, separator into easier to manage location
write.table(resp,file=paste(project,"_fc.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

print(sum(bad1 & bad2))

# TODO: are these bads for the 2 genes or the 2 replicates?  Former ok, latter not
good<-!bad1 & !bad2

rge<-range(fc)
# TODO: where do the 0.005s in the break seq come from?  Do they need to be set in params?
hh<-hist(fc[!allbad],breaks=seq(rge[1]-0.005,rge[2]+0.005,by=0.005),main="construct fitness",col="grey80",border=FALSE,xlab=expression(f["c"]))
d<-density(fc[!allbad],bw=0.01)
lines(d$x,d$y*sum(hh$counts)*0.005,col="black")

rug(fc[!allbad])
abline(h=0)

plot_scatterplots()

resp<-data.frame(probes,p_rank)
# TODO: refactor file suffix, separator into easier to manage location
write.table(resp,file=paste(project,"_p.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

resy<-data.frame(genes,f,f_sd)
# TODO: refactor col names, file suffix, separator into easier to manage location
colnames(resy)<-c("gene","f","sd")
write.table(resy[-1,],file=paste(project,"_f.txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)

# TODO: since this is always on diagonal, should we refactor to be rug plot?

rge<-range(f)*1.1 #with some margin
plot(f[-1],f[-1],pch=16,cex=0.6,col="blue",xlim=rge,ylim=rge,xlab=expression(f["g"]),ylab=expression(f["g"]),main="single-gene fitness")
abline(v=0,h=0,col="#000066")
# TODO: Do I need to worry about these 2s?
text(f[-1],f[-1],cex=0.6,labels=genes[-1],pos=((rank(f[-1]) %% 2)+1)*2,offset=0.2)