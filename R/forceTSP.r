
forceTSPNext <- function(distance_matix, path, distance)
{
	pointNum <- dim(distance_matix)[1]
	pathLength <- length(path)
	minDistance <- -1
	bestPath <- c()
	#message("pointnum = ", pointNum, " path = ", path, ", length = ", pathLength, ", distance = ", distance)
	if (pointNum == pathLength)
	{
		#message("return")
		return(c(distance+distance_matix[path[pathLength], path[1]], path))
	}
	tmp <- c()
	for (i in 1: pointNum)
	{
		if (length(which(path == i)) == 0){
			tmp <- forceTSPNext(distance_matix, c(path, i), distance+distance_matix[path[pathLength], i])
			#message("return tmp length is", length(tmp))
			if (minDistance < 0 | minDistance > tmp[1])
			{
				minDistance <- tmp[1]
				bestPath <- tmp[2: length(tmp)]
			}
		}
	}
	#message("return")
	return(c(minDistance, bestPath))
}

forceTSP <- function(distance_matix)
{
	pointNum <- dim(distance_matix)
	tmp <- c()
	tmp <- forceTSPNext(distance_matix, c(1), 0)
	minDistance <- tmp[1]
	bestPath <- tmp[2: length(tmp)]
	bestPath <- c(bestPath, 1)
	return(c(minDistance, bestPath))
}
