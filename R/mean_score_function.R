# get mean score

averageScoreall <- function(data, ccgenes)
{

  rname <- row.names(data)
  
  yy <- match(tolower(ccgenes[,1]), tolower(rname))
  
  xx <- which(is.na(yy))
  
  if (length(xx) != 0){
    ccgenes <- ccgenes[-xx, ]
  }
  
  itsindex <- which(!is.na(match(tolower(rname), tolower(ccgenes[,1]))))
  
  tmp <- data[itsindex, ]
  
  #average <- rowMeans(tmp)
  
  #tmp <- apply(tmp, 2, function(x) x - average)
  
  G1index <- which(ccgenes[,2] == "G1")
  G1Sindex = which(ccgenes[,2] == "G1/S")
  Sindex <- which(ccgenes[,2] == "S")
  G2index <- which(ccgenes[,2] == "G2")
  G2Mindex <- which(ccgenes[,2] == "G2/M") 
  Mindex = which(ccgenes[,2] == "M")
  
  xx <- match(tolower(row.names(tmp)), tolower(ccgenes[,1]))
  
  yy <- which(!is.na(match(xx, G1index)))
  zz <- which(!is.na(match(xx, Sindex)))
  uu <- which(!is.na(match(xx, G2Mindex)))
  vv <- which(!is.na(match(xx, G1Sindex)))
  ww <- which(!is.na(match(xx, G2index)))
  tt <- which(!is.na(match(xx, Mindex)))
  
  G1Score <- apply(tmp, 2, function(x) mean(x[yy]))
  SScore <- apply(tmp, 2, function(x) mean(x[zz]))
  G2MScore <- apply(tmp, 2, function(x) mean(x[uu]))
  G1SScore <- apply(tmp, 2, function(x) mean(x[vv]))
  G2Score <- apply(tmp, 2, function(x) mean(x[ww]))
  MScore <- apply(tmp, 2, function(x) mean(x[tt]))
  
  result.average <- data.frame(G1Score = G1Score, G1SScore = G1SScore, SScore = SScore, G2Score = G2Score, G2MScore = G2MScore, MScore = MScore)
  
  return(result.average)
}