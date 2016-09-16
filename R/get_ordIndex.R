# get the order
# 
# you can choose thread you want to run the order (threadnum)


source("bestEnsembleComplexTSP.r")

get_ordIndex <- function(test_exp, threadnum, cls_num)
{
  result <- bestEnsembleComplexTSP(test_exp = test_exp, TSPFold = 2, beginNum = 10, endNum = cls_num, threads = threadnum)
  ord <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], ]
  ordIndex <- order(ord)
  return(ordIndex)
}
