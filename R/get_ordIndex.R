source("bestEnsembleComplexTSP.r")

get_ordIndex <- function(test_exp, threadnum)
{
  result <- bestEnsembleComplexTSP(test_exp = test_exp, TSPFold = 2, beginNum = 10, endNum = dim(test_exp)[1], threads = threadnum)
  ord <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], ]
  ordIndex <- order(ord)
  return(ordIndex)
}
