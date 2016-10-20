source("get_hmm.R")
library(doParallel)

get_start <- function(bayes_score, mean_score, ordIndex, cls_num, rdata = NULL, nthread = 3)
{
	cl <- makeCluster(nthread)
	registerDoParallel(cl)
	le = length(ordIndex)
	start = 0
	p = -Inf
	for(i in 1:le)
	{
		myord = c (i:1, le:(i+1))
		if (i == le)
			myord = c(le:1)
		
	    re = get_hmm_order(bayes_score, mean_score, ordIndex, cls_num, myord, rdata)
		if (max(re$p) > p)
		{
			start = i
			p = max(re$p)
		}    	
	}

	#p = apply(log_lk, 2, max)
	#start = which(p == max(p))

	return(start)
}
