normalDis <- function(val)
{
  len = dim(val)[1]
  t = dim(val)[2]
  m = apply(val, 2, mean)
  v = apply(val, 2, function(x) var(x)*(len-1)/len)
  result <- data.frame(m = m, v = sqrt(v))
  return(result)
}

calcEmissionProbs <- function(ndresult, nowval, M)
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

mulP <- function(plist)
{
  result = 1
  for (i in 1:length(plist))
  {
    result = result * plist[i]
  }
  
  return(result)
}


forwardPro <- function(transProb, nd, mypi, obval, M)
{
  pmat = matrix(nrow = dim(obval)[1], ncol = M)
  c = rep(0, dim(obval)[1])
  nowval = obval[1,]
  nowemissionp = calcEmissionProbs(ndresult = nd, nowval = nowval, M = M)
  sum = 0
  for(i in 1:M)
  {
    pmat[1,i] = mulP(nowemissionp[i,])*mypi[i]
    sum = sum + pmat[1,i]
  }
  c[1] = 1 / sum
  for(i in 1:M)
  {
    pmat[1,i] = c[1] * pmat[1,i]
  }
  
  
  for (i in 2:dim(obval)[1])
  {
    nowval = obval[i,]
    nowemissionp = calcEmissionProbs(ndresult = nd, nowval = nowval, M = M)
    sum = 0
    for(j in 1:M)
    {
      calcp = mulP(nowemissionp[j,])
      x = 0
      for (k in 1:M){x = x + pmat[i-1,k]*transProb[k,j]}
      pmat[i, j] = x*calcp
      sum = sum + pmat[i, j]
    }
    c[i] = 1 / sum
    for(j in 1:M)
    {
      pmat[i,j] = c[i] * pmat[i,j]
    }
  }
  
  return(list(pmat = pmat, c = c))
}


backwardPro <- function(transProb, nd, mypi, obval, M, c = c)
{
  pmat = matrix(nrow = dim(obval)[1], ncol = M)
  for(i in 1:M)
  {
    pmat[dim(obval)[1],i] = c[dim(obval)[1]]
  }
  
  for (i in dim(obval)[1]:2)
  {
    nowval = obval[i,]
    nowemissionp = calcEmissionProbs(ndresult = nd, nowval = nowval, M = M)
    calcp = apply(nowemissionp, 1, mulP)
    for(j in 1:M)
    {
      x = 0
      for (k in 1:M){x = x + calcp[k]*transProb[j,k]*pmat[i,k]}
      pmat[i-1, j] = x*c[i-1]
    }
  }
  
  return(pmat)
}

gammaPro <- function(fPro, bPro, transProb, nd, mypi, obval, M, c)
{
  mr = array(0, dim = c(M, M, dim(obval)[1]))
  r = matrix(0, dim(obval)[1], M)
  for (t in 1:(dim(obval)[1]-1))
  {
    nowval = obval[t+1,]
    nowemissionp = calcEmissionProbs(ndresult = nd, nowval = nowval, M = M)
    
    for (i in 1:M)
    {
      for (j in 1:M)
      {
        mr[i,j,t] = fPro[t,i]*transProb[i,j]*mulP(nowemissionp[j,])*bPro[t+1,j]
      }
      r[t,i] = fPro[t,i]*bPro[t,i] / c[t]
    }
  }
  
  t = dim(obval)[1]
  for (i in 1:M)
  {
    r[t,i] = fPro[t, i] * bPro[t,i] / c[t]
  }
  
  return(list(mr = mr, r = r))
}


calclogP <- function(c)
{
  p = 0
  for (i in 1:length(c))
  {
    p = p + log(c[i])
  }
  return(-p)
}


myBW <- function(transProb, nd, mypi, obval, N, M)
{
  q = rep(0, N)
  for(k in 1:N)
  {
    message(k)
    fPro = forwardPro(transProb = transProb, nd = nd, mypi = mypi, obval = obval, M = M)
    bPro = backwardPro(transProb = transProb, nd = nd, mypi = mypi, obval = obval, M = M, c = fPro$c)
    
    gammaP = gammaPro(fPro = fPro$pmat, bPro = bPro, transProb = transProb, nd = nd, mypi = mypi, obval = obval, M = M, c = fPro$c)
    
    
    newmypi = rep(0, M)
    for (i in 1:M)
    {
      newmypi[i] = gammaP$r[1,i]
    }
    
    
    newtransPro = matrix(0, M, M)
    for (i in 1:M)
    {
      for (j in 1:M)
      {
        numer = 0
        denom = 0
        for (t in 1:(dim(obval)[1]-1))
        {
          numer = numer + gammaP$mr[i,j,t]
          denom = denom + gammaP$r[t,i]
        }
        newtransPro[i,j] = numer/denom
      }
    }
   
    
    u = rep(0, dim(nd)[1])
    v = rep(0, dim(nd)[1])
    for (i in 1:M)
    {
      for (j in 1:dim(obval)[2])
      {
        numer = 0
        denom = 0
        for (t in 1:dim(obval)[1])
        {
          numer = numer + gammaP$r[t,i]*obval[t,j]
          denom = denom + gammaP$r[t,i]
        }
        u[(i-1)*dim(obval)[2] + j] = numer/denom
        numer = 0
        denom = 0
        for (t in 1:dim(obval)[1])
        {
          numer = numer + gammaP$r[t,i]*(obval[t,j]-u[(i-1)*dim(obval)[2] + j])^2
          denom = denom + gammaP$r[t,i]
        }
        v[(i-1)*dim(obval)[2] + j] = numer/denom
        if (v[(i-1)*dim(obval)[2] + j] < 0.0001){v[(i-1)*dim(obval)[2] + j] = 0.0001}
      }
    }
    newnd = cbind(u, sqrt(v))
    
    
    transProb = newtransPro
    nd = newnd
    mypi = newmypi
    q[k] = calclogP(fPro$c)
  }
  
  return(list(transProb = transProb, nd = nd, mypi = mypi, q = q))
}