source("HMM.R")
source("baum_welch_scale.R")

get_bw_three <- function(bayes_score, mean_score, G1.id, S.id, G2M.id, ordIndex, fob)
{
  transPro = matrix(0, 3, 3)
  transPro[1,1] = 0.5
  transPro[1,2] = 0.5
  transPro[2,2] = 0.5
  transPro[2,3] = 0.5
  transPro[3,3] = 0.5
  transPro[3,1] = 0.5
  
  mypi = c(1/3, 1/3, 1/3, 1/3)
  
  pred2 = bayes_score[ordIndex, ]
  result2 = mean_score[ordIndex, ]
 
  G1val = data.frame(x1 = result2$G1Score[G1.id], x2 = result2$G1SScore[G1.id],x3 = result2$SScore[G1.id],
                     x4 = result2$G2Score[G1.id], x5 = result2$G2MScore[G1.id], x6 = result2$MScore[G1.id],
                     x7 = pred2$G1.score[G1.id], x8 = pred2$S.score[G1.id], x9 = pred2$G2M.score[G1.id])
  
  Sval = data.frame(x1 = result2$G1Score[S.id], x2 = result2$G1SScore[S.id],x3 = result2$SScore[S.id],
                    x4 = result2$G2Score[S.id], x5 = result2$G2MScore[S.id], x6 = result2$MScore[S.id],
                    x7 = pred2$G1.score[S.id], x8 = pred2$S.score[S.id], x9 = pred2$G2M.score[S.id])
  
  G2Mval = data.frame(x1 = result2$G1Score[G2M.id], x2 = result2$G1SScore[G2M.id],x3 = result2$SScore[G2M.id],
                      x4 = result2$G2Score[G2M.id], x5 = result2$G2MScore[G2M.id], x6 = result2$MScore[G2M.id],
                      x7 = pred2$G1.score[G2M.id], x8 = pred2$S.score[G2M.id], x9 = pred2$G2M.score[G2M.id])
  
  G1nd = normalDis(G1val)
  Snd = normalDis(Sval)
  G2Mnd = normalDis(G2Mval)
  ndresult <- rbind(G1nd, Snd, G2Mnd)
  
  if (fob == 1)
  {
    obval = data.frame(x1 = (result2$G1Score),x2 = (result2$G1SScore),x3 = (result2$SScore),
                       x4 = (result2$G2Score), x5 = (result2$G2MScore), x6 = (result2$MScore),
                       x7 = (pred2$G1.score), x8 = (pred2$S.score), x9 = (pred2$G2M.score))
  }
  else
  {
    obval = data.frame(x1 = rev(result2$G1Score),x2 = rev(result2$G1SScore),x3 = rev(result2$SScore),
                       x4 = rev(result2$G2Score), x5 = rev(result2$G2MScore), x6 = rev(result2$MScore),
                       x7 = rev(pred2$G1.score), x8 = rev(pred2$S.score), x9 = rev(pred2$G2M.score))
  }
  
  obval = as.matrix(obval)
  
  myresult = myBW(transProb = transPro, nd = ndresult, mypi = mypi, obval = obval, 20, M = 3)
  rs1 = myviterbi(obval = obval, transProb = myresult$transProb, ndresult = myresult$nd, mypi = myresult$mypi, M = 3)
  
  nowr = order(rs1$rmat[length(ordIndex),])[3]
  rr = c(nowr)
  for (i in (length(ordIndex)-1):1)
  {
    rr = c(rr,rs1$pmat[i,nowr])
    nowr = rs1$pmat[i,nowr]
  }
  rr <- rev(rr)
  rr[1] = order(rs1$rmat[1,])[3]
  
  return(bw_order = rr)
}


get_bw_four <- function(bayes_score, mean_score, G0.id, G1.id, S.id, G2M.id, ordIndex, fob)
{
  transPro = matrix(0, 4, 4)
  transPro[1,1] = 0.5
  transPro[1,2] = 0.5
  transPro[2,2] = 0.5
  transPro[2,3] = 0.5
  transPro[3,3] = 0.5
  transPro[3,4] = 0.5
  transPro[4,4] = 0.5
  transPro[4,1] = 0.5
  
  mypi = c(1/4, 1/4, 1/4, 1/4)
  
  pred2 = bayes_score[ordIndex, ]
  result2 = mean_score[ordIndex, ]
  
  G0val = data.frame(x1 = result2$G1Score[G0.id], x2 = result2$G1SScore[G0.id],x3 = result2$SScore[G0.id],
                     x4 = result2$G2Score[G0.id], x5 = result2$G2MScore[G0.id], x6 = result2$MScore[G0.id],
                     x7 = pred2$G1.score[G0.id], x8 = pred2$S.score[G0.id], x9 = pred2$G2M.score[G0.id])
  
  G1val = data.frame(x1 = result2$G1Score[G1.id], x2 = result2$G1SScore[G1.id],x3 = result2$SScore[G1.id],
                     x4 = result2$G2Score[G1.id], x5 = result2$G2MScore[G1.id], x6 = result2$MScore[G1.id],
                     x7 = pred2$G1.score[G1.id], x8 = pred2$S.score[G1.id], x9 = pred2$G2M.score[G1.id])
  
  Sval = data.frame(x1 = result2$G1Score[S.id], x2 = result2$G1SScore[S.id],x3 = result2$SScore[S.id],
                    x4 = result2$G2Score[S.id], x5 = result2$G2MScore[S.id], x6 = result2$MScore[S.id],
                    x7 = pred2$G1.score[S.id], x8 = pred2$S.score[S.id], x9 = pred2$G2M.score[S.id])
  
  G2Mval = data.frame(x1 = result2$G1Score[G2M.id], x2 = result2$G1SScore[G2M.id],x3 = result2$SScore[G2M.id],
                      x4 = result2$G2Score[G2M.id], x5 = result2$G2MScore[G2M.id], x6 = result2$MScore[G2M.id],
                      x7 = pred2$G1.score[G2M.id], x8 = pred2$S.score[G2M.id], x9 = pred2$G2M.score[G2M.id])
  
  G0nd = normalDis(G0val)
  G1nd = normalDis(G1val)
  Snd = normalDis(Sval)
  G2Mnd = normalDis(G2Mval)
  ndresult <- rbind(G0nd, G1nd, Snd, G2Mnd)
  
  if (fob == 1)
  {
    obval = data.frame(x1 = (result2$G1Score),x2 = (result2$G1SScore),x3 = (result2$SScore),
                       x4 = (result2$G2Score), x5 = (result2$G2MScore), x6 = (result2$MScore),
                       x7 = (pred2$G1.score), x8 = (pred2$S.score), x9 = (pred2$G2M.score))
  }
  else
  {
    obval = data.frame(x1 = rev(result2$G1Score),x2 = rev(result2$G1SScore),x3 = rev(result2$SScore),
                       x4 = rev(result2$G2Score), x5 = rev(result2$G2MScore), x6 = rev(result2$MScore),
                       x7 = rev(pred2$G1.score), x8 = rev(pred2$S.score), x9 = rev(pred2$G2M.score))
  }
  
  obval = as.matrix(obval)
  
  myresult = myBW(transProb = transPro, nd = ndresult, mypi = mypi, obval = obval, 20, M = 4)
  rs1 = myviterbi(obval = obval, transProb = myresult$transProb, ndresult = myresult$nd, mypi = myresult$mypi, M = 4)
  
  nowr = order(rs1$rmat[length(ordIndex),])[4]
  rr = c(nowr)
  for (i in (length(ordIndex)-1):1)
  {
    rr = c(rr,rs1$pmat[i,nowr])
    nowr = rs1$pmat[i,nowr]
  }
  rr <- rev(rr)
  rr[1] = order(rs1$rmat[1,])[4]
  
  return(bw_result = rr)
}