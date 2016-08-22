# there are funcions to do HMM

transProb <- function(statlist)
{
  smat = matrix(0,nrow = 3, ncol = 3)
  x = statlist[1]
  for (i in 2:length(statlist)){
    y = statlist[i]
    smat[x,y] = smat[x,y] + 1
    x = y
  }
  
  pmat = matrix(0,nrow = 3, ncol = 3)
  pmat[1,1] = smat[1,1] / (smat[1,1] + smat[1,2])
  pmat[1,2] = smat[1,2] / (smat[1,1] + smat[1,2])
  pmat[2,2] = smat[2,2] / (smat[2,2] + smat[2,3])
  pmat[2,3] = smat[2,3] / (smat[2,2] + smat[2,3])
  pmat[3,3] = smat[3,3] / (smat[3,1] + smat[3,3])
  pmat[3,1] = smat[3,1] / (smat[3,1] + smat[3,3])
  
  return(pmat)
}

normalDis <- function(val)
{
  len = dim(val)[1]
  t = dim(val)[2]
  m = apply(val, 2, mean)
  v = apply(val, 2, function(x) var(x)*(len-1)/len)
  result <- data.frame(m = m, v = sqrt(v))
  return(result)
}


calcEmissionProbs <- function(ndresult, nowval, M = M)
{
  result = matrix(nrow = M, ncol = length(nowval))
  for (i in 0:(M-1))
  {
    for (j in 1:length(nowval))
    {
      result[i+1,j] = dnorm(nowval[j], ndresult[i*length(nowval)+j,1], ndresult[i*length(nowval)+j,2])
    }
  }
  
  return(result)
}


vmulP <- function(plist)
{
  result = 0
  for (i in 1:length(plist))
  {
    result = result + log2(plist[i])
  }
  
  return(result)
}

myviterbi <- function(obval, transProb, ndresult, mypi, M)
{
  rmat = matrix(nrow = dim(obval)[1], ncol = M)
  pmat = matrix(nrow = dim(obval)[1], ncol = M)
  nowval = obval[1,]
  nowemissionp = calcEmissionProbs(ndresult = ndresult, nowval = nowval, M = M)
  for(i in 1:M)
  {
    rmat[1,i] = vmulP(nowemissionp[i,])+log2(mypi[i]) 
    pmat[1,i] = 0
  }
  
  for (i in 2:dim(obval)[1])
  {
    nowval = obval[i,]
    nowemissionp = calcEmissionProbs(ndresult = ndresult, nowval = nowval, M = M)
    for(j in 1:M)
    {
      calcp = vmulP(nowemissionp[j,])
      x = c()
      for (k in 1:M){x = c(x, rmat[i-1,k]+log2(transProb[k,j]))}
      pre = max(x)
      rmat[i,j] = calcp + pre
      pmat[i,j] = which(x == pre)[1]
    }
  }
  
  return(list(rmat = rmat, pmat = pmat))
}