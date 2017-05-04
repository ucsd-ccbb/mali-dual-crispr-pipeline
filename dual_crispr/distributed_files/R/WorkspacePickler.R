getLsVarMask <- function(env, maskFunction) {
  result = sapply(ls(envir = env, all.names = TRUE), function(x) maskFunction(x))
  return(result)
}


getRelevantNames <- function(env, excludedNames = NULL){
  nonFunctionMask = getLsVarMask(env, function(x) !is.function(get(x)))
  notExcludedMask = rep(TRUE, length(nonFunctionMask))
  if (!is.null(excludedNames)) {
    notExcludedMask = getLsVarMask(env, function(x) !(x %in% excludedNames))
  }

  combinedMask = nonFunctionMask & notExcludedMask
  relevantVarNames = ls(envir = env, all.names = TRUE)[combinedMask]
  return(relevantVarNames)
}


idAndCompareVariable <- function(varName, currVarValue, oldEnv, oldVarNames) {
  outputCurrObj = TRUE # default assumption
  inOldEnv = FALSE

  if (!is.null(oldEnv)) {
    if (varName %in% oldVarNames) {
      inOldEnv = TRUE
      oldObjValue <- get(varName, envir = oldEnv)
      if (identical(oldObjValue, currVarValue)) {
        outputCurrObj = FALSE
      } else {
        print(sprintf('CHANGED  %s', varName) )
      }
    }
  }

  if (!inOldEnv) print(sprintf('NEW %s', varName) )

  if (outputCurrObj) {
    coercible = TRUE
    tryCatch({
      as.data.frame(currVarValue)
    }, error = function(e) {
      coercible <<- FALSE
    })
    if (!coercible) {
      print(sprintf('%s not coercible ', varName) )
      outputCurrObj = FALSE
    }
  }

  return(outputCurrObj)
}


saveWorkspaceVariables <- function(outputDir, newTimeptId, oldTimeptId=NULL, excludedVarNames = NULL) {
  getTimeptFpPrefix <- function(timeptId) {
    return(file.path(outputDir, timeptId))
  }

  getTimeptRdataFp <- function(timeptId){
    timeptFpPrefix = getTimeptFpPrefix(timeptId)
    return(paste0(timeptFpPrefix,".RData"))
  }

  save.image(file = getTimeptRdataFp(newTimeptId))
  newTimeptFpPrefix =  getTimeptFpPrefix(newTimeptId)
  newTimeptEnv <- globalenv()
  objs <- getRelevantNames(newTimeptEnv, excludedVarNames)

  oldTimeptEnv = NULL
  oldObjs = NULL
  if (!is.null(oldTimeptId)) {
    oldTimeptEnv = new.env()
    load(getTimeptRdataFp(oldTimeptId), envir = oldTimeptEnv)
    oldObjs = getRelevantNames(oldTimeptEnv, excludedVarNames)
  }

  for (obj in objs) {
    currObjValue <- get(obj, envir = newTimeptEnv)
    writeOut = idAndCompareVariable(obj, currObjValue, oldTimeptEnv, oldObjs)
    if (writeOut) write.csv(currObjValue, file = paste0(newTimeptFpPrefix, "_", obj, '.csv'))
  }
}
