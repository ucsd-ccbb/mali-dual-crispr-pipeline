# prepData=function(input_filename, nt){
#   X<-read.table(input_filename,sep="\t",header=TRUE)
#
#   #preliminary preparations of the input data frame
#   data<-data.matrix(X[,6:(5+2*nt)])
#   good<-(X$geneA != X$geneB) #reject any constructs with two 0's
#   goodX<-X[good,] #the 0-0 constructs are gone
#   nn<-sum(good) #this many constructs
#
#   cpA<-as.character(goodX$probeA)
#   ix<-grep("NonTargeting",cpA)
#   cpA[ix]<-paste("0",cpA[ix],sep="") #this puts NonTargeting probes at the beginning of alphabetically sorted order
#
#   cpB<-as.character(goodX$probeB)
#   ix<-grep("NonTargeting",cpB)
#   cpB[ix]<-paste("0",cpB[ix],sep="")
#
#   pswitch<-cpA>cpB #need to switch?
#   phold<-cpA[pswitch]
#   cpA[pswitch]<-cpB[pswitch]
#   cpB[pswitch]<-phold #cpA and cpB are always in alphabetical order, cpA < cpB
#   probes<-sort(unique(c(cpA,cpB))) #entire probe set in alphabetical order
#   return(probes)
#
#   nprobes<-length(probes)
#
#   cgA<-as.character(goodX$geneA)
#   cgB<-as.character(goodX$geneB)
#   genes<-sort(unique(cgA)) #should be 74 "genes"
#   n<-length(genes) # n = 74 if doing it by genes or 222 if doing it by probe
#   mm<-n*(n-1)/2
#
#   gswitch<-cgA>cgB #need to switch?
#   ghold<-cgA[gswitch]
#   cgA[gswitch]<-cgB[gswitch]
#   cgB[gswitch]<-ghold
#
#   gA_gB<-paste(cgA,cgB,sep="_")
#   pA_pB<-paste(cpA,cpB,sep="_")
#   goodX<-data.frame(goodX,cgA,cgB,gA_gB) #now gA_gB is ordered so that gA < gB
#
#   gooddata<-data.matrix(goodX[,6:(5+2*nt)])
#   gooddata[gooddata==0]<-1 #pseudocounts
#   abundance<-apply(gooddata,2,sum)
#   y<-t(log2(t(gooddata)/abundance)) #log2 frequencies
#   #end data prep
#   return(list(nn=nn, probes=probes, nprobes=nprobes, genes=genes, n=n, pA_pB=pA_pB, y=y))
# }
# 
# getMeasurementColsAsMatrix = function(countsDf, numTimepoints) {
#   data <- data.matrix(countsDf[, 6:(5 + 2 * numTimepoints)])
#   return(data)
# }
# 
# excludeConstructsForSameTwoGenes = function(countsDf) {
#   good <-
#     (countsDf$geneA != countsDf$geneB) #reject any constructs with two 0's
#   goodX <- countsDf[good, ] #the 0-0 constructs are gone
#   return(goodX)
# }
# 
# renameNontargetingContentsForCol = function(filteredCountsDf, colName) {
#   cpA <- as.character(filteredCountsDf[[colName]])
#   ix <- grep("NonTargeting", cpA)
#   cpA[ix] <-
#     paste("0", cpA[ix], sep = "") #this puts NonTargeting probes at the beginning of alphabetically sorted order
#   return(cpA)
# }
# 
# orderItemsInVectorPairsAlphabetically = function(firstVector, secondVector) {
#   pswitch <- firstVector > secondVector #need to switch?
#   phold <- firstVector[pswitch]
#   firstVector[pswitch] <- secondVector[pswitch]
#   secondVector[pswitch] <-
#     phold #firstVector and secondVector are always in alphabetical order, firstVector < secondVector
#   return(list(vector1 = firstVector, vector2 = secondVector))
# }
# 
# getSortedUniqueItems = function(firstVector, secondVector) {
#   sortedUniques <-
#     sort(unique(c(firstVector, secondVector))) #entire probe set in alphabetical order
#   return(sortedUniques)
# }
# 
# prepDataChanged = function(input_filename, nt) {
#   countsDf <- read.table(input_filename, sep = "\t", header = TRUE)
#   goodX = excludeConstructsForSameTwoGenes(countsDf)
#   nn = nrow(goodX)
#   
#   cpA = renameNontargetingContentsForCol(goodX, "probeA")
#   cpB = renameNontargetingContentsForCol(goodX, "probeB")
#   probeVectors = orderItemsInVectorPairsAlphabetically(cpA, cpB)
#   probes <-
#     getSortedUniqueItems(probeVectors$vector1, probeVectors$vector2)
#   return(probes)
#   
#   cgA = renameNontargetingContentsForCol(goodX, "geneA")
#   cgB = renameNontargetingContentsForCol(goodX, "geneB")
#   geneVectors = orderItemsInVectorPairsAlphabetically(cgA, cgB)
#   genes <-
#     getSortedUniqueItems(geneVectors$vector1, geneVectors$vector2)
#   
#   gA_gB <- paste(cgA, cgB, sep = "_")
#   pA_pB <- paste(cpA, cpB, sep = "_")
#   goodX <-
#     data.frame(goodX, cgA, cgB, gA_gB) #now gA_gB is ordered so that gA < gB
#   
#   gooddata <- getMeasurementColsAsMatrix(goodX, numTimepoints)
#   gooddata[gooddata == 0] <- 1 #pseudocounts
#   abundance <- apply(gooddata, 2, sum)
#   y <- t(log2(t(gooddata) / abundance)) #log2 frequencies
#   #end data prep
#   
#   return(list(nn = nn, probes = probes, nprobes = nprobes, genes = genes, n = n, pA_pB = pA_pB, y = y))
# }
