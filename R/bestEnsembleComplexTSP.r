# load(file = "~/Florian/highCorFindFlo.RData")
# message("Task Begin! open ~/Florian/highCorFindFlo.RData")
# beginNum = 10
# endNum = 198
# output = TRUE
# debug = TRUE
# name = ""
# threads = 1
source("highCorFind.r")

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

bestEnsembleComplexTSP <- function(test_exp, stageIdxSmp = NULL, TSPFold = 2, beginNum = 10, endNum = 198, output = TRUE, debug = FALSE, name = "", threads = 1, clustMethod = "GMM", base_cycle_range = c(6:9), step_size = 2, max_loop = NULL, max_num = 300)
{
	if (debug == TRUE)
		message("Task Begin! ")
	require(doParallel)
	cl <- makeCluster(threads)
	
	stageTest <- !is.null(stageIdxSmp)
	
	if (stageTest)
	{
		G1IdxSmp <- stageIdxSmp[["G1"]]
		SIdxSmp <- stageIdxSmp[["S"]]
		G2MIdxSmp <- stageIdxSmp[["G2M"]]
		stageToNum <- rep(c(1,2,3),c(length(G1IdxSmp),length(SIdxSmp),length(G2MIdxSmp)))
		cellNum <- length(G1IdxSmp) + length(SIdxSmp) + length(G2MIdxSmp)
	}

	# initialize
  #resultLst <- matrix(0, nrow = (endNum - beginNum + 1), ncol = nrow(test_exp))
	tmpResultLst <- matrix(0, 4, ncol = nrow(test_exp))
	#normalizedResultLst <- matrix(0, nrow = (endNum - beginNum + 1), ncol = nrow(test_exp))
	#ensembleResultLst <- matrix(0, nrow = (endNum - beginNum + 1), ncol = nrow(test_exp))
	
	#print(tmpResultLst)
	beginNum <- max(base_cycle_range)
	if(endNum > 500){
	 endNum <- max_num
	}
	
  if(is.null(max_loop))
	{
    s1 <- (endNum - beginNum) %% step_size
	  total_rows <- round2(((endNum - beginNum - s1) / step_size), 0) + 1 #+ s1
	  resultLst <- matrix(0, nrow = total_rows, ncol = nrow(test_exp))
	  normalizedResultLst <- matrix(0, nrow = total_rows, ncol = nrow(test_exp))
	  ensembleResultLst <- matrix(0, nrow = total_rows, ncol = nrow(test_exp))
	  seq_list <- seq(from = beginNum+step_size, to = endNum, by = step_size)
	}
	else if (!is.null(max_loop))
	{
	  s1 <- (endNum - max_loop) %% step_size
	  max_single_step <- max_loop - beginNum
	  total_rows <- (max_loop - beginNum) + round2(((endNum - max_loop - s1) / step_size), 0) + 1 # + s1 
	  resultLst <- matrix(0, nrow = total_rows, ncol = nrow(test_exp))
	  normalizedResultLst <- matrix(0, nrow = total_rows, ncol = nrow(test_exp))
	  ensembleResultLst <- matrix(0, nrow = total_rows, ncol = nrow(test_exp))
	  seq_list <- c((beginNum+1):max_loop, seq(from = max_loop+step_size, to = endNum, by = step_size))
	}


	
	#print(total_rows)
	#print(dim(resultLst))
	
	highCorLst <- c()
	finalHighCorLst <- c()
  #pause(0.1)
	# caculate the base score
  base_cycle_range <-base_cycle_range
	for (clust_num_EM in base_cycle_range) #in 7: 10)
	{
		i <- clust_num_EM - (min(base_cycle_range)-1)#6
		highCorTour <- highCorFind(clust_num_EM,test_exp,stageIdxSmp,TSPFold = TSPFold, TSPMethod = "force", clustMethod = clustMethod, debug = debug)
		if (debug == TRUE)
		message(clust_num_EM, " Clusts, high Cor is ", highCorTour$highCor)
		tmpResultLst[i, ] <- highCorTour$result
	}
  #print(tmpResultLst)
  #pause(0.1)
	tmpCorLst <- c()
	ri <- tmpResultLst[1, ]
	first_j <- 2*base_cycle_range[1] + 1
	second_k <- 2*base_cycle_range[2] + 1
	third_l <- 2*base_cycle_range[3] + 1
	
	j1 <- base_cycle_range[2]
	k1 <- base_cycle_range[3]
	l1 <- base_cycle_range[4]
	
	for (j in 0:first_j)#0: 15)
	{
	  #print(c("first loop now is : ", j))
	  if (j < j1)#8)
	  {
	    #rj <- (tmpResultLst[2, ] -0.5 + j) %% 8
	    rj <- rotateResult(j1, tmpResultLst[2, ], rotate = j)#8
	  }
	  else
	  {
	    #rj <- (j - tmpResultLst[2, ] - 0.5) %% 8
	    rj <- rotateResult(j1, tmpResultLst[2, ], rotate = j-j1, back = TRUE)
	  }
	  for (k in 0: second_k)#17)
	  {
	    #print(c("second loop now is : ", k))
	    if (k < k1)#9)
	    {
	      #rk <- (tmpResultLst[3, ] -0.5 + k) %% 9
	      rk <- rotateResult(k1, tmpResultLst[3, ], rotate = k)#9, tmpResultLst[3, ], rotate = k)
	      
	    }
	    else
	    {
	      #rk <- (k - tmpResultLst[3, ] - 0.5) %% 9
	      rk <- rotateResult(k1, tmpResultLst[3, ], rotate = k-k1, back = TRUE)#9, tmpResultLst[3, ], rotate = k-9, back = TRUE)
	    }
	    for (l in 0: third_l)#19)
	    {
	      #print(c("third loop now is : ", l))
	      if (l < l1)#10)
	      {
	        #rl <- (tmpResultLst[4, ] - 0.5 + l) %% 10 
	        rl <- rotateResult(l1, tmpResultLst[4, ], rotate = l)#10, tmpResultLst[4, ], rotate = l)
	      }
	      else
	      {
	        #rl <- (l - tmpResultLst[4, ] - 0.5) %% 10
	        rl <- rotateResult(l1, tmpResultLst[4, ], rotate = l-l1, back = TRUE)#10, tmpResultLst[4, ], rotate = l-10, back = TRUE)
	      }
	      tmpCor <- cor(ri, rj)
	      tmpCor <- tmpCor + cor(ri, rk)
	      tmpCor <- tmpCor + cor(ri, rl)
	      tmpCor <- tmpCor + cor(rj, rk)
	      tmpCor <- tmpCor + cor(rj, rl)
	      tmpCor <- tmpCor + cor(rl, rk)
	      tmpCorLst <- c(tmpCorLst, tmpCor)
	    }
	  }
	}

	highNum <- which(tmpCorLst == max(tmpCorLst))[1]-1
	
	#print(highNum)
	
	j2 <- 2*base_cycle_range[2]
	k2 <- 2*base_cycle_range[3]
	l2 <- 2*base_cycle_range[4]
	
	j3 <- base_cycle_range[1]+1
	k3 <- base_cycle_range[2]+1
	l3 <- base_cycle_range[3]+1
	
	l <- highNum %% l2#20 
	highNum <- (highNum - l) / l2#20 
	k <- highNum %% k2#18
	highNum <- (highNum - k) / k2#18
	j <- highNum %% j2#16

	#print(l)
	#print(k)
	#print(j)
	#print(max(tmpCorLst))
	tmpResultLst[1, ] <- rotateResult(j1-1, tmpResultLst[1, ], normalize = 1)#7, tmpResultLst[1, ], normalize = 1)
	if (j < j3)#8)
	{
		tmpResultLst[2, ] <- rotateResult(j3, tmpResultLst[2, ], rotate = j, normalize = 1)#8, tmpResultLst[2, ], rotate = j, normalize = 1)
	}else
	{
		tmpResultLst[2, ] <- rotateResult(j3, tmpResultLst[2, ], rotate = j - j3, normalize = 1, back = TRUE)#8, tmpResultLst[2, ], rotate = j - 8, normalize = 1, back = TRUE)
	}
	if (k < k3)#9)
	{
		tmpResultLst[3, ] <- rotateResult(k3, tmpResultLst[3, ], rotate = k, normalize = 1)#9, tmpResultLst[3, ], rotate = k, normalize = 1)
	}else
	{
		tmpResultLst[3, ] <- rotateResult(k3, tmpResultLst[3, ], rotate = k - k3, normalize = 1, back = TRUE)#9, tmpResultLst[3, ], rotate = k - 9, normalize = 1, back = TRUE)
	}
	if (l < 10)#10)
	{
		tmpResultLst[4, ] <- rotateResult(10, tmpResultLst[4, ], rotate = l, normalize = 1)#10, tmpResultLst[4, ], rotate = l, normalize = 1)
	}else
	{
		tmpResultLst[4, ] <- rotateResult(10, tmpResultLst[4, ], rotate = l - 10, normalize = 1, back = TRUE)#10, tmpResultLst[4, ], rotate = l - 10, normalize = 1, back = TRUE)
	}
	baseResult <- colMeans(tmpResultLst)
	
  #print(baseResult)
  
	normalizedResultLst[1, ] <- baseResult
	ensembleResultLst[1, ] <- baseResult
	if (stageTest)
	{
		finalHighCorLst <- c(cor(baseResult, stageToNum))
		highCorLst <- finalHighCorLst
		if (debug == TRUE)
			message("mean Cor of 7-10 is ", finalHighCorLst)
	}
	#pause(0.1)
	

	if (threads <= 1)
	{
		# sequential processing
		resultLst[1, ] <- colMeans(tmpResultLst)
		#for (clust_num_EM in (beginNum + 1): endNum)
		#x1 = 2
		i <- 2
		for (clust_num_EM in seq_list)#seq(from = beginNum+step_size, to = endNum, by = step_size))
		{
			#i <- clust_num_EM - beginNum + 1
		  #i <- 2# - beginNum + 1
		  
			highCorTour <- highCorFind(clust_num_EM, test_exp, stageIdxSmp, TSPFold = TSPFold, clustMethod = clustMethod, debug = debug)
			if (debug == TRUE && stageTest)
				message(clust_num_EM, " Clusts, high Cor is ", highCorTour$highCor)
			
			resultLst[i, ] <- highCorTour$result
			i <- i + 1
			#x1 = x1 + 1
		}
	}
	else
	{
		# parallel processing
	  print("make cluster")
		cl <- makeCluster(threads)
		print("register cluster")
		registerDoParallel(cl)
		print("finished")
		#resultLst <- foreach(clust_num_EM = (beginNum + 1): endNum, .combine='rbind') %dopar% 
		resultLst <- foreach(clust_num_EM = seq_list, .combine = 'rbind') %dopar% #seq(from = beginNum+step_size, to = endNum, by = step_size), .combine = 'rbind') %dopar%
		{
			source("highCorFind.r")
			highCorTour <- highCorFind(clust_num_EM, test_exp, stageIdxSmp, TSPFold = TSPFold, clustMethod = clustMethod, debug = debug)
			return(highCorTour$result)
		}
		stopCluster(cl)
		closeAllConnections()
		resultLst <- rbind(baseResult, resultLst)
		if (debug == TRUE)
			message(Sys.time())
	}
	assign("res_new", resultLst, .GlobalEnv)
	
	
	#pause(0.1)
	# rotate the time score, to make the corelation with the base score be the best
	#for (clust_num_EM in (beginNum + 1): endNum)
	#y1 = 0
	
	i <- 2
	for (clust_num_EM in seq_list)#seq(from = beginNum+step_size, to = endNum, by = step_size))
	{
		#i <- clust_num_EM - beginNum + 1
		#i <- 2#clust_num_EM - beginNum + step_size
		if (stageTest)
		{
			highCorLst <- c(highCorLst, cor(resultLst[i, ], stageToNum))
		}
		resulti <- resultLst[i, ]
		moveCorLst <- c()
		backMoveCorLst <- c()
		#tmpCor <- c()
		#backTmpCor <- c()

		for (k in 0: (clust_num_EM - 1))
		{
			# tmp <- (resulti + k - 0.5) %% clust_num_EM
			tmp <- rotateResult(clust_num_EM, resulti, rotate = k)
			moveCorLst <- c(moveCorLst, cor(ensembleResultLst[1, ], tmp))
			#tmpCor <- c(tmpCor, cor(stageToNum, tmp))
			# tmp <- (clust_num_EM - resulti + 1 + k - 0.5) %% clust_num_EM
			tmp <- rotateResult(clust_num_EM, resulti, rotate = k, back = TRUE)
			backMoveCorLst <- c(backMoveCorLst, cor(ensembleResultLst[1, ], tmp))
			#backTmpCor <- c(backTmpCor, cor(stageToNum, tmp))
		}
		moveMaxID <- which(moveCorLst == max(moveCorLst))[1] - 1
		backMoveMaxID <- which(backMoveCorLst == max(backMoveCorLst))[1] - 1
		if (max(moveCorLst) >= max(backMoveCorLst))
		{
			# normalizedResultLst[i, ] <- ((resulti + moveMaxID - 1 - 0.5) %% clust_num_EM) / clust_num_EM 
			normalizedResultLst[i, ] <- rotateResult(clust_num_EM, resulti, rotate = moveMaxID, normalize = 1)
			if (debug == TRUE)
		    message("moved ", moveMaxID, ", Cor is ", max(moveCorLst))
			
		}
		else
		{
			# normalizedResultLst[i, ] <- ((clust_num_EM - resulti + backMoveMaxID - 0.5) %% clust_num_EM) / clust_num_EM
			normalizedResultLst[i, ] <- rotateResult(clust_num_EM, resulti, rotate = backMoveMaxID, back = TRUE, normalize = 1)
			if (debug == TRUE)
				message("back moved ", backMoveMaxID, ", Cor is ", max(backMoveCorLst))
			
		}
		#print(i)
		i <- i + 1
	}
	if (debug == TRUE)
		message(Sys.time())
	
	#for (clust_num_EM in (beginNum + 1): endNum)
	#print(dim(normalizedResultLst))
	#print(normalizedResultLst[1,])
	assign("norm_new", normalizedResultLst, .GlobalEnv)
	i <- 2
	for (clust_num_EM in seq_list)#seq(from = beginNum+step_size, to = endNum, by = step_size))
	{
	  #i <- 2#clust_num_EM - beginNum + 1
	  #print(clust_num_EM)
	  meanResults <- colMeans(matrix(complex(argument = (normalizedResultLst[1: i, ] * 2 * pi)), nrow = i))
	  meanResults <- Arg(meanResults) / 2 / pi
	  tmp <-which(meanResults < 0)
	  meanResults[tmp] <- meanResults[tmp] + 1
	  ensembleResultLst[i, ] <- meanResults
	  if (stageTest)
	  {
	    sortedResult <- sort(ensembleResultLst[i, ])
	    tmpHighCorLst <- c()
	    for (j in 1: cellNum)
	    {
	      tmpResult <- ensembleResultLst[i, ] - sortedResult[j]
	      tmp <- which(tmpResult < 0)
	      tmpResult[tmp] <- tmpResult[tmp] + 1
	      tmpHighCorLst <- c(tmpHighCorLst, cor(tmpResult, stageToNum))
	      tmpResult <- 1 - (ensembleResultLst[i, ] - sortedResult[j])
	      tmp <- which(tmpResult < 0)
	      tmpResult[tmp] <- tmpResult[tmp] + 1
	      tmpHighCorLst <- c(tmpHighCorLst, cor(tmpResult, stageToNum))
	    }
	    finalHighCor <- max(tmpHighCorLst)
	    if (debug == TRUE)
	      message("final ", i, ": moved ", which(tmpHighCorLst == finalHighCor) / 2 - 1)
	    #finalHighCor <- cor(meanResults, stageToNum)
	    if (debug == TRUE)
	      message("mean HighCor is ", finalHighCor)
	    finalHighCorLst <- c(finalHighCorLst, finalHighCor)
	  }
	  i <- i + 1
	  #print("second: ")
	  #print(i)
	}
	

	if (debug == TRUE)
		message("print in ", paste("bestEnsembleComplexTSP", beginNum, "-", endNum, name, ".pdf"))
		
	if (output == TRUE)
	{
		save(finalHighCorLst, highCorLst, resultLst, normalizedResultLst, ensembleResultLst, file = paste("bestEnsembleComplexTSP", beginNum, "-", endNum, name, ".RData"))
		
		pdf(paste("bestEnsembleComplexTSP", beginNum, "-", endNum, name, ".pdf"))
		plot(beginNum: endNum, finalHighCorLst, xlab = "Clust Num", ylab = "highCor", type = "o", ylim = c(0, 1))
		lines(beginNum: endNum, highCorLst, col = "red")
		dev.off()
	}
	
	return(list(finalHighCorLst = finalHighCorLst, ensembleResultLst = ensembleResultLst))
}
