##What is reCAT?

reCAT is a modelling framework for unsynchronized single-cell transcriptome data that can reconstruct a high-resolution cell cycle time-series. It is mainly based on traveling salesman problem (TSP), Gaussian mixture model (GMM) and hidden Markov model (HMM). We developed an R based software package which is easy to use. The performance is relatively accurate and stable. Thanks for using!

Software by Zehua Liu, Huazhe Lou and Hao Wang. reCAT is explained in more detail in the accompanying publication. 

##Philosophy

Two fundamental assumptions underlie the reCAT approach: that 1) different cell cycle phases form a cycle, and 2) the change of transcriptome profile from one phase to the next should be monotonic, increasing with time span widths. 

##Installation


##How to use

reCAT is easy to use. Now we use ola_mES_2i in /data as example

####1.preparatory work

when you use our tools, you should install some packages first, the package list is as follows:

***ggplot2***

***doParallel***

***mclust***

***cluster***

***TSP***

####2.input data

in reCAT, there are some  requirements for the input data. 

1. the genes must be in ***cyclebaseGeneList***
2. 

*ola_mES_2i.RData* in */data* is an example

####3.get order

when you preprocessing your test data, you can get its order(cell's time series) easily with ***get_ordIndex*** function. In this function, there are two parameters, one is the input data, the other is thread number, so maybe you can choose a large thread number like 20 to speed it up. 
for example:

	source("get_ordIndex.R")
	load("../data/ola_mES_2i.RData")
	ordIndex <- get_ordIndex(test_exp, 10)

####4.get bayes-score and mean score
in reCAT, there are two scores to 
for example:

	source("get_score.R")
	score_result <- getScore(t(test_exp))
you can use the following two orders to get scores

	score_result$bayes_score
	score_result$mean_score

####5.plot1
plot with order you get:

	source("plot.R")
	plot_bayes(score_result$bayes_score, ordIndex)
	plot_mean(score_result$mean_score, ordIndex)
the result is like follows:

<div align="center">
<img src="./pic/ola_2i_bayes.png" width = "300" height = "200" alt="ola_2i_bayes"/>
<img src="./pic/ola_2i_mean.png" width = "300" height = "200" alt="ola_2i_mean"/>
</div>

####6.HMM
for example:
	
	source("get_bw.R")
	load("../data/ola_mES_2i_ordIndex.RData")
	load("../data/ola_mES_2i_region.RData")
	hmm_result <- get_bw_three(score_result$bayes_score, score_result$mean_score, ordIndex, cls_num = 3, fob = 0)

####7.plot2
plot with HMM result:
	
	source("plot.R")
	load("../data/ola_mES_2i_hmm.RData")
	plot_bayes(score_result$bayes_score, ordIndex, hmm_result, hmm_order, 1)
	plot_mean(score_result$mean_score, ordIndex, hmm_result, hmm_order, 1)
the result is like follows:

<div align="center">
<img src="./pic/ola_2i_bayes_hmm.png" width = "300" height = "200" alt="ola_2i_bayes_hmm"/>
<img src="./pic/ola_2i_mean_hmm.png" width = "300" height = "200" alt="ola_2i_mean_hmm"/>
</div>