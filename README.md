##How to use
----------
reCAT is easy to use. Now we use ola_mES_2i in /data as example

####input data

in reCAT, there are some  requirements for the input data. 

1. the genes must be in ***cyclebaseGeneList***
2. 

*ola_mES_2i.RData* in */data* is an example

####get order

when you preprocessing your test data, you can get its order(cell's time series) easily with ***get_ordIndex*** function. In this function, there are two parameters, one is the input data, the other is thread number, so maybe you can choose a large thread number to speed it up like 20. 
for example:

	source("get_ordIndex.R")
	load("../data/ola_mES_2i.RData")
	ordIndex <- get_ordIndex(test_exp, 10)

####get bayes-score and mean score
in reCAT, there are two scores to 
for example:

	source("get_score.R")
	score_result <- getScore(t(test_exp))
you can use the following two orders to get scores

	score_result$bayes_score
	score_result$mean_score

####HMM
for example:
	
	source("get_bw.R")
	load("../data/ola_mES_2i_ordIndex.RData")
	load("../data/ola_mES_2i_region.RData")
	hmm_result <- get_bw_three(score_result$bayes_score, score_result$mean_score, G1.id, S.id, G2M.id, ordIndex, 0)

####plot
plot with order you get:

	source("plot.R")
	plot_bayes_base(score_result$bayes_score, ordIndex)
	plot_mean_base(score_result$mean_score, ordIndex)
plot with HMM result:
	
	source("plot.R")
	load("../data/ola_mES_2i_hmm.RData")
	plot_bayes_bw(score_result$bayes_score, ordIndex, hmm_result, hmm_order)
	plot_mean_bw(score_result$mean_score, ordIndex, hmm_result, hmm_order)