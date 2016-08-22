 
# load(file = "~/Florian/highCorFindFlo.RData")
# message("Task Begin! open ~/Florian/highCorFindFlo.RData")
# beginNum = 10
# endNum = 198
# output = TRUE
# debug = TRUE
# name = ""
# threads = 1
	
source("highCorFind.r")
bestEnsembleComplexTSP <- function(test_exp, stageIdxSmp = NULL, TSPFold = 2, beginNum = 10, endNum = 198, output = TRUE, debug = FALSE, name = "", threads = 1, clustMethod = "GMM")
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
	resultLst <- matrix(0, nrow = (endNum - beginNum + 1), ncol = nrow(test_exp))
	tmpResultLst <- matrix(0, 4, ncol = nrow(test_exp))
	normalizedResultLst <- matrix(0, nrow = (endNum - beginNum + 1), ncol = nrow(test_exp))
	ensembleResultLst <- matrix(0, nrow = (endNum - beginNum + 1), ncol = nrow(test_exp))
	highCorLst <- c()
	finalHighCorLst <- c()

	# caculate the base score
	for (clust_num_EM in 7: 10)
	{
		i <- clust_num_EM - 6
		highCorTour <- highCorFind(clust_num_EM,test_exp,stageIdxSmp,TSPFold = TSPFold, TSPMethod = "force", clustMethod = clustMethod, debug = debug)
		if (debug == TRUE)
		message(clust_num_EM, " Clusts, high Cor is ", highCorTour$highCor)
		tmpResultLst[i, ] <- highCorTour$result
	}

	tmpCorLst <- c()
	ri <- tmpResultLst[1, ]
	for (j in 0: 15)
	{
		if (j < 8)
		{
			#rj <- (tmpResultLst[2, ] -0.5 + j) %% 8
			rj <- rotateResult(8, tmpResultLst[2, ], rotate = j)
		}
		else
		{
			#rj <- (j - tmpResultLst[2, ] - 0.5) %% 8
			rj <- rotateResult(8, tmpResultLst[2, ], rotate = j-8, back = TRUE)
		}
		for (k in 0: 17)
		{
			if (k < 9)
			{
				#rk <- (tmpResultLst[3, ] -0.5 + k) %% 9
				rk <- rotateResult(9, tmpResultLst[3, ], rotate = k)
				
			}
			else
			{
				#rk <- (k - tmpResultLst[3, ] - 0.5) %% 9
				rk <- rotateResult(9, tmpResultLst[3, ], rotate = k-9, back = TRUE)
			}
			for (l in 0: 19)
			{
				if (l < 10)
				{
					#rl <- (tmpResultLst[4, ] - 0.5 + l) %% 10 
					rl <- rotateResult(10, tmpResultLst[4, ], rotate = l)
				}
				else
				{
					#rl <- (l - tmpResultLst[4, ] - 0.5) %% 10
					rl <- rotateResult(10, tmpResultLst[4, ], rotate = l-10, back = TRUE)
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
	l <- highNum %% 20 
	highNum <- (highNum - l) / 20 
	k <- highNum %% 18
	highNum <- (highNum - k) / 18
	j <- highNum %% 16

	tmpResultLst[1, ] <- rotateResult(7, tmpResultLst[1, ], normalize = 1)
	if (j < 8)
	{
		tmpResultLst[2, ] <- rotateResult(8, tmpResultLst[2, ], rotate = j, normalize = 1)
	}else
	{
		tmpResultLst[2, ] <- rotateResult(8, tmpResultLst[2, ], rotate = j - 8, normalize = 1, back = TRUE)
	}
	if (k < 9)
	{
		tmpResultLst[3, ] <- rotateResult(9, tmpResultLst[3, ], rotate = k, normalize = 1)
	}else
	{
		tmpResultLst[3, ] <- rotateResult(9, tmpResultLst[3, ], rotate = k - 9, normalize = 1, back = TRUE)
	}
	if (l < 10)
	{
		tmpResultLst[4, ] <- rotateResult(10, tmpResultLst[4, ], rotate = l, normalize = 1)
	}else
	{
		tmpResultLst[4, ] <- rotateResult(10, tmpResultLst[4, ], rotate = l - 10, normalize = 1, back = TRUE)
	}
	baseResult <- colMeans(tmpResultLst)

	normalizedResultLst[1, ] <- baseResult
	ensembleResultLst[1, ] <- baseResult
	if (stageTest)
	{
		finalHighCorLst <- c(cor(baseResult, stageToNum))
		highCorLst <- finalHighCorLst
		if (debug == TRUE)
			message("mean Cor of 7-10 is ", finalHighCorLst)
	}

	if (threads <= 1)
	{
		# sequential processing
		resultLst[1, ] <- colMeans(tmpResultLst)
		for (clust_num_EM in (beginNum + 1): endNum)
		{
			i <- clust_num_EM - beginNum + 1
			highCorTour <- highCorFind(clust_num_EM, test_exp, stageIdxSmp, TSPFold = TSPFold, clustMethod = clustMethod, debug = debug)
			if (debug == TRUE && stageTest)
				message(clust_num_EM, " Clusts, high Cor is ", highCorTour$highCor)
			
			resultLst[i, ] <- highCorTour$result
		}
	}
	else
	{
		# parallel processing
		cl <- makeCluster(threads)
		registerDoParallel(cl)
		resultLst <- foreach(clust_num_EM = (beginNum + 1): endNum, .combine='rbind') %dopar%
		{
			source("highCorFind.r")
			highCorTour <- highCorFind(clust_num_EM, test_exp, stageIdxSmp, TSPFold = TSPFold, clustMethod = clustMethod, debug = debug)
			return(highCorTour$result)
		}
		stopCluster(cl)
		resultLst <- rbind(baseResult, resultLst)
		if (debug == TRUE)
			message(Sys.time())
	}
	# retote the time score, to make the corelation with the base score be the best
	for (clust_num_EM in (beginNum + 1): endNum)
	{
		i <- clust_num_EM - beginNum + 1
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
	}
	if (debug == TRUE)
		message(Sys.time())
	
	for (clust_num_EM in (beginNum + 1): endNum)
	{
		i <- clust_num_EM - beginNum + 1
		
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