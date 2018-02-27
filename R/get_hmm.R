# there are two funcitons to get HMM, main function is get_HMM_order
# 
# the parameters are:
#   
#   bayes_score : calculate before
#   mean_score : calculate before
#   cls_num : class number you want
#   rdata : include (G0), G1, S, G2M region. 
#           when you get order, you can plot use plot_bayes_base or plot_mean_base function in plot.R
#           so you can judge where is G1, where is S and where is G2/M, the parameter is the region you judge
#           rdata = t(data.frame(c(2,5), c(10,15), c(20,25))), means from 2 to 5 is G1, 10 to 15 is S and 20 to 25 is G2M                  
#   ordIndex : calculate before
#   fob : when the cells are ordered in ordIndex, the time series is forward(=1) or backward(=0) 

# when we get HMM, we should do baum welch first to get the origrin transition probability and emission probability
# 
                      


source("HMM.R")
source("baum_welch_scale_log.R")

get_region <- function(list, len)
{
  if (list[1] < list[2])
  {
    return(c(list[1]:list[2]))
  }
  else
  {
    return(c(list[1]:len, 1:list[2]))
  }
}

get_hmm_order <- function(bayes_score, mean_score, ordIndex, cls_num, myord, rdata = NULL)
{
  pred2 = bayes_score[ordIndex, ]
  result2 = mean_score[ordIndex, ]
  
  if (cls_num == 3)
  {
    transPro = matrix(0, 3, 3)
    transPro[1,1] = 0.5
    transPro[1,2] = 0.5
    transPro[2,2] = 0.5
    transPro[2,3] = 0.5
    transPro[3,3] = 1
    #transPro[3,1] = 0
    
    mypi = c(1/3, 1/3, 1/3)
    
    if(is.null(rdata))
    {
      rdata = matrix(nrow = 3, ncol = 2)
      print("please input G1's region like 2 5")
      rdata[1,] = scan()
      print("please input S's region like 2 5")
      rdata[2,] = scan()
      print("please input G2M's region like 2 5")
      rdata[3,] = scan()
    }
    
    G1.id = get_region(rdata[1,], length(ordIndex))
    S.id = get_region(rdata[2,], length(ordIndex))
    G2M.id = get_region(rdata[3,], length(ordIndex))
    
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
  }
  
  if (cls_num == 4)
  {
    transPro = matrix(0, 4, 4)
    transPro[1,1] = 0.5
    transPro[1,2] = 0.5
    transPro[2,2] = 0.5
    transPro[2,3] = 0.5
    transPro[3,3] = 0.5
    transPro[3,4] = 0.5
    transPro[4,4] = 1
    #transPro[4,1] = 0.5
    
    mypi = c(1/4, 1/4, 1/4, 1/4)
    
    if(is.null(rdata))
    {
      rdata = matrix(nrow = 4, ncol = 2)
      print("please input G0's region like 2 5")
      rdata[1,] = scan()
      print("please input G1's region like 2 5")
      rdata[2,] = scan()
      print("please input S's region like 2 5")
      rdata[3,] = scan()
      print("please input G2M's region like 2 5")
      rdata[4,] = scan()
    }
    
    G0.id = get_region(rdata[1,], ordIndex)
    G1.id = get_region(rdata[2,], ordIndex)
    S.id = get_region(rdata[3,], ordIndex)
    G2M.id = get_region(rdata[4,], ordIndex)
    
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
  }
  
 
  obval = data.frame(x1 = (result2$G1Score)[myord],x2 = (result2$G1SScore)[myord],x3 = (result2$SScore)[myord],
                     x4 = (result2$G2Score)[myord], x5 = (result2$G2MScore)[myord], x6 = (result2$MScore)[myord],
                     x7 = (pred2$G1.score)[myord], x8 = (pred2$S.score)[myord], x9 = (pred2$G2M.score)[myord])
  
  obval = as.matrix(obval)

  #print("good!")
  #save(transPro,ndresult,mypi,obval,M=cls_num,file="C:/experiment/BW_test.RData")  
  iter_max = 30
  myresult = myBW(A = transPro, nd_para = ndresult, mypi = mypi, ob_value = obval, iter_max = iter_max)
   
  rs1 = myviterbi(obval = obval, transProb = myresult$transProb, ndresult = myresult$nd, mypi = myresult$mypi, M = cls_num)
    
  nowr = order(rs1$rmat[length(ordIndex),])[cls_num]
  rr = c(nowr)
  for (i in (length(ordIndex)-1):1)
  {
    rr = c(rr,rs1$pmat[i,nowr])
    nowr = rs1$pmat[i,nowr]
  }
  rr <- rev(rr)
  rr[1] = order(rs1$rmat[1,])[cls_num]

  return(rr)
}