#01. Compute eexp(x)
eexp <- function(x)
{
  if(is.nan(x)==TRUE)
  {
    return(0)
  }
  else
  {
    return(exp(x))
  }
}

#02. Compute eln(x)
eln <- function(x)
{
  if(x==0)
  {
    return(NaN)
  }
  else if (x > 0)
  {
    return(log(x))
  }
  else
  {
    stop('negative imput error')
  }
}

#03. Compute elnsum(eln(x),eln(y)) ==  ln(x+y)
#For (x+y)
elnsum <- function(eln_x, eln_y)
{
  if(is.nan(eln_x)==TRUE | is.nan(eln_y)==TRUE)
  {	
	if(is.nan(eln_x)==TRUE)
    {
      return(eln_y)
    }
	else
    {
      return(eln_x)
	}
  }
  else
  {
    if(eln_x > eln_y)
	{
	  return(eln_x+eln(1+exp(eln_y-eln_x)))
	}
	else
	{
	  return(eln_y+eln(1+exp(eln_x-eln_y)))
	}
  }
}

#04. Compute elnproduct(eln(x),eln(y)) == ln(x*y)
#For (x*y)
elnproduct <- function(eln_x, eln_y)
{
  if(is.nan(eln_x)==TRUE | is.nan(eln_y)==TRUE)
  {
    return(NaN)
  }
  else
  {
    return(eln_x+eln_y)
  }
}

#05. Compute eln(αt(i)) for all states Si and observations Ot
eln_Alpha <- function(eln_A, eln_pi, eln_b)
{
  M <- dim(eln_A)[1]
  N <- dim(eln_b)[1]
  eln_alpha = matrix(0, nrow = N , ncol = M)
  for(i in 1:M)
  {
    eln_alpha[1,i] <- elnproduct(eln_pi[i], eln_b[1,i])
  }
  for(t in 2:N)
  {
    for(j in 1:M)
	{
	  logalpha <- NaN
	  for(i in 1:M)
	  {
	    logalpha <- elnsum(logalpha, elnproduct(eln_alpha[t-1,i],eln_A[i,j]))
	  }
	  eln_alpha[t,j] <- elnproduct(logalpha, eln_b[t,j])
	}
  }
  return(eln_alpha)
}

#06. Compute eln(βt(i)) for all states Si and observations Ot
eln_Beta <- function(eln_A, eln_b)
{
  M <- dim(eln_A)[1]
  N <- dim(eln_b)[1]
  eln_beta = matrix(0, nrow = N, ncol = M)
  for(i in 1:M)
  {
    eln_beta[N,i] <- 0
  }
  for(t in (N-1):1)
  {
    for(i in 1:M)
	{
	  logbeta <- NaN
	  for(j in 1:M)
	  {
	    logbeta <- elnsum(logbeta, elnproduct(eln_A[i,j], elnproduct(eln_b[t+1,j], eln_beta[t+1,j])))
	  }
	  eln_beta[t,i] <- logbeta
	}
  }
  return(eln_beta)
}

#07. Compute eln(γt(i)) for all states Si and observations Ot
eln_Gamma <- function(eln_alpha,eln_beta)
{
  N <- dim(eln_alpha)[1]
  M <- dim(eln_alpha)[2]
  eln_gamma <- matrix(0, N, M)
  for(t in 1:N)
  {
    normalizer <- NaN
	for(i in 1:M)
	{
	  eln_gamma[t,i] <- elnproduct(eln_alpha[t,i],eln_beta[t,i])
	  normalizer <- elnsum(normalizer,eln_gamma[t,i])
	}
	for(i in 1:M)
	{
	  eln_gamma[t,i] <- elnproduct(eln_gamma[t,i],-normalizer)
	}
  }
  return(eln_gamma)  
}

#08. Compute eln(ξt(i, j)) for all state pairs Si and Sj and observations Ot
eln_Xi <- function(eln_A, eln_b, eln_alpha, eln_beta)
{
  N <- dim(eln_alpha)[1]
  M <- dim(eln_alpha)[2]
  eln_xi <- array(0, dim = c(N, M, M))
  for(t in 1:(N-1))
  {
    normalizer <- NaN
	for(i in 1:M)
	{
	  for(j in 1:M)
	  {
	    eln_xi[t,i,j] <- elnproduct(eln_alpha[t,i],elnproduct(eln_A[i,j],elnproduct(eln_b[t+1,j], eln_beta[t+1,j])))
	    normalizer <- elnsum(normalizer,eln_xi[t,i,j])
	  }
	}
	for(i in 1:M)
	{
	  for(j in 1:M)
	  {
	    eln_xi[t,i,j] <- elnproduct(eln_xi[t,i,j],-normalizer)
	  }
	}
  }
  return(eln_xi)
}

#09. Compute π(i), the estimated probability of starting in state Si
get_Pi <- function(eln_gamma)
{
  M <- dim(eln_gamma)[2]
  pi_new <- rep(0,M)
  for(i in 1:M)
  {
    pi_new[i] <- eexp(eln_gamma[1,i])
  }
  return(pi_new)
}

#10. Compute a(i,j) , the estimated probability of transitioning from state Si to Sj
get_A <- function(eln_gamma, eln_xi)
{
  N <- dim(eln_gamma)[1]
  M <- dim(eln_gamma)[2]
  A_new <- matrix(0, M, M)
  for(i in 1:M)
  {
    for(j in 1:M)
    {
      numerator <- NaN
      denominator <- NaN
      for(t in 1:(N-1))
      {
        numerator <- elnsum(numerator, eln_xi[t,i,j])
	    denominator <- elnsum(denominator, eln_gamma[t,i])
      }
      A_new[i,j] <- eexp(elnproduct(numerator,-denominator))
	}
  }
  return(A_new)
}

#11. Compute parameters for the emission probability
get_Ndpara <- function(eln_gamma, ob_value)
{
  ob_dim <- dim(ob_value)[2]
  N <- dim(eln_gamma)[1]
  M <- dim(eln_gamma)[2]
  mu <- rep(0, M*ob_dim)
  va <- rep(0, M*ob_dim)
  for (i in 1:M)
  {
    for (j in 1:ob_dim)
    {
      numerator = NaN
      denominator = NaN
      for (t in 1:N)
      {
	    #abs() for the observations
        numerator <- elnsum(numerator, elnproduct(eln_gamma[t,i],eln(abs(ob_value[t,j]))))
        denominator <- elnsum(denominator, eln_gamma[t,i])
      }
      mu[(i-1)*ob_dim + j] <- sign(ob_value[t,j])* eexp(elnproduct(numerator,-denominator))
      
      numerator = NaN
      denominator = NaN
      for (t in 1:N)
      {
        numerator <- elnsum(numerator, elnproduct(eln_gamma[t,i],eln((ob_value[t,j]-mu[(i-1)*ob_dim + j])^2)))
        denominator <- elnsum(denominator, eln_gamma[t,i])
      }
      va[(i-1)*ob_dim + j] <- eexp(elnproduct(numerator,-denominator))
      #print(va[(i-1)*ob_dim + j])
	  #考虑要不要这一段
      if (va[(i-1)*ob_dim + j] < 0.0001)
      {
	    va[(i-1)*ob_dim + j] = 0.0001
      }
    }
  }
  nd_new = cbind(mu, sqrt(va))
  return(nd_new)
}

#12. Compute the current emission probability
calc_lnb <- function(nd_para, ob_value)
{
  N <- dim(ob_value)[1]
  ob_dim <- dim(ob_value)[2]
  M <- dim(nd_para)[1]/ob_dim
  eln_b = matrix(0, nrow = N, ncol = M)
  for(t in 1:N)
  {
    for(i in 1:M)
	{
	  lnprob <- 0
	  for (j in 1:ob_dim)
      {
        lnprob <- elnproduct(lnprob, eln(dnorm(ob_value[t,j], nd_para[(i-1)*ob_dim+j,1], nd_para[(i-1)*ob_dim+j,2])))
	  }
	  eln_b[t,i] <- lnprob
	}	
  }
  return(eln_b)
}

#13. Main Baum-Welch funciton
myBW <- function(A, nd_para, mypi, ob_value, iter_max)
{
  N <- dim(ob_value)[1]
  ob_dim <- dim(ob_value)[2]
  M <- dim(nd_para)[1]/ob_dim
  q = rep(-Inf, iter_max)
  A_new <- A
  nd_new <- nd_para
  pi_new <- mypi
  
  for(k in 1:iter_max)
  {
  
    #Do logarithm preparation
    eln_A <- matrix(0, M, M)
    eln_pi <- rep(0, M)
    for(i in 1:M)
      for(j in 1:M)
        eln_A[i,j] <- eln(A_new[i,j])
    for(i in 1:M)
      eln_pi[i] <- eln(pi_new[i])
    eln_b <- calc_lnb(nd_new, ob_value)
    #print(eln_b)

	#Formal Baum-Welch
	eln_alpha <- eln_Alpha(eln_A, eln_pi, eln_b)
    eln_beta <- eln_Beta(eln_A, eln_b)
    eln_gamma <- eln_Gamma(eln_alpha,eln_beta)
	eln_xi <- eln_Xi(eln_A, eln_b, eln_alpha, eln_beta)
	pi_new <- get_Pi(eln_gamma)
	#print(pi_new)
	#print(eln_beta)
	#print(eln_gamma)
	#print(eln_xi)
    A_new <- get_A(eln_gamma, eln_xi)
	#print(A_new)
    nd_new <- get_Ndpara(eln_gamma, ob_value)
    #print(nd_new)
    
	fP_obv = NaN
	for(l in 1:M)
	{
	  fP_obv <- elnsum(fP_obv, eln_alpha[N,l])
	}
	#print(fP_obv)
    q[k] <- fP_obv
	
	if(k!=1)
	  if(q[k]-q[k-1] < 1e-6)
	    break
  }
  
  return(list(transProb = A_new, nd = nd_new, mypi = pi_new, q = q))
}







