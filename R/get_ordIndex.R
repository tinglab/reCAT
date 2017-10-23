# get the order
# 
# you can choose thread you want to run the order (threadnum)


source("bestEnsembleComplexTSP.r")


get_ordIndex <- function(test_exp, threadnum, step_size = 2, base_cycle_range = c(6:9), max_loop = NULL)
{
  result <- bestEnsembleComplexTSP(test_exp = test_exp, TSPFold = 2, beginNum = 10, endNum = dim(test_exp)[1], threads = threadnum, base_cycle_range = base_cycle_range, step_size = step_size, max_loop = max_loop, max_num = 300)
  ord <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], ]
  #print(ord)
  ordIndex <- order(ord)
  return(ordIndex)
}


#get_ordIndex_original <- function(test_exp, threadnum)
#{
#  result <- bestEnsembleComplexTSP(test_exp = test_exp, TSPFold = 2, beginNum = 10, endNum = dim(test_exp)[1], threads = threadnum)
#  ord <- result$ensembleResultLst[dim(result$ensembleResultLst)[1], ]
#  ordIndex <- order(ord)
#  return(ordIndex)
#}
