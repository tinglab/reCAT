source("bestEnsembleComplexTSP.r")

get_ordIndex <- function(test_exp)
{
  result <- bestEnsembleComplexTSP(test_exp = test_exp, TSPFold = 2, beginNum = 10, endNum = dim(test_exp)[1], threads = 20)
  ord <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], ]
  ordIndex <- order(ord)
  return(ordIndex)
}
