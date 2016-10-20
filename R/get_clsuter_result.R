source("get_score.R")

get_cluster_result <- function(test_exp, ensembleResultLst, resultLst, cls_num)
{
	if (cls_num <= 10)
	{
		print("the cluster number is too small!")
		return(0)
	}
	else
	{
		rownum = cls_num-10+1

		ordIndex2 = ensembleResultLst[rownum, ]
		EM_result = resultLst[rownum, ]
		cls_mean = c()
		for(i in c(1:cls_num))
		{
		  id = which(EM_result == i)
		  tmp = sum(complex(argument = ordIndex2[id] * 2 * pi - pi))
		  cls_mean = cbind(cls_mean, tmp)	
		}
		meanResults <- Arg(cls_mean) / 2 / pi
		ordIndex3 = order(meanResults)

		cls_mean = c()
		for(i in c(1:cls_num))
		{
		  id = which(EM_result == i)
		  
		  if (length(id) == 1)
		  {
		    tmp = test_exp[id, ]
		  }
		  else
		  {
		    tmp = apply(test_exp[id, ], 2, mean)
		  }
		  cls_mean = cbind(cls_mean, tmp)
		}

		score_r = get_score(cls_mean)

		return(list(cls = EM_result, bayes_score = score_r$bayes_score, mean_score = score_r$mean_score, ordIndex = ordIndex3))							
	}
}