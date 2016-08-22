# there are functions that do Naive bayes
# 
# the result is bayes_score


getPairs <- function(data, G1.marker.pairs, 
                     S.marker.pairs, G2M.marker.pairs, genes.training){
  genes.test<-row.names(data)
  
  markers<-unique(as.character(unlist(c((G1.marker.pairs[,1:2]),
                                        (S.marker.pairs[,1:2]),
                                        (G2M.marker.pairs[,1:2])))))
  
  if(all(markers%in%genes.test)==T){#checking if markers are in the list of genes of the testing dataset
    
    genes.list<-intersect(genes.training, genes.test)
    G1.markers<-G1.marker.pairs
    S.markers<-S.marker.pairs
    G2M.markers<-G2M.marker.pairs
    
  }else{
    
    genes.list<-intersect(genes.training, genes.test)
    genes.to.remove<-markers[!markers%in%genes.test]    #remove markers that are not in the dataset
    id<-t(apply(G1.marker.pairs[,1:2],1,function(x) match(genes.to.remove,x)))    
    id<-apply(id, 1, function(x) all(is.na(x)) )
    id<-which(id==F)
    if(length(id)!=0){
      G1.markers<-G1.marker.pairs[-id,]
    }else{
      G1.markers<-G1.marker.pairs
    }
    
    id<-t(apply(S.marker.pairs[,1:2],1,function(x) match(genes.to.remove,x)))    
    id<-apply(id, 1, function(x) all(is.na(x)) )
    id<-which(id==F)
    if(length(id)!=0){
      S.markers<-S.marker.pairs[-id,]
    }else{
      S.markers<-S.marker.pairs
    }
    
    id<-t(apply(G2M.marker.pairs[,1:2],1,function(x) match(genes.to.remove,x)))    
    id<-apply(id, 1, function(x) all(is.na(x)) )
    id<-which(id==F)
    if(length(id)!=0){
      G2M.markers<-G2M.marker.pairs[-id,]
    }else{
      G2M.markers<-G2M.marker.pairs
    }    
    
  }
  
  print(sprintf("Number of G1 pairs: %d", nrow(G1.markers)))
  print(sprintf("Number of S pairs: %d", nrow(S.markers)))
  print(sprintf("Number of G2M pairs: %d", nrow(G2M.markers)))
  
  markers <- rbind(G1.markers, S.markers, G2M.markers)
  
  return(markers)
}

vector.single <- function(cell, markers)
{
  test <- unlist(cell[markers[,1]] - cell[markers[,2]])
  
  v <- rep(-1, length(test))
  
  index <- which(test > 0)
  
  v[index] <- 1
  
  return(v)
}

getVector <- function(data, markers)
{
  markers[,1]<-match(as.character(markers[,1]),row.names(data))
  markers[,2]<-match(as.character(markers[,2]),row.names(data))
  
  result <- apply(data, 2, function(x) vector.single(x, markers))
  
  return(result)  
}

cellscore <- function(cell, p, q)
{
  result <- 0
  for (i in 1:length(cell)){
    if (cell[i] == 1)
    {
      result <- result + logb(p[i,1])
    }else{
      result <- result + logb(p[i,2])
    }
  }
  
  result <- result + logb(q)
  
  return(result)
}

Naive <- function(data.training, G1.l, S.l, G2M.l, len, data.test)
{
  P.G1 <- G1.l/len
  P.S <- S.l/len
  P.G2M <- G2M.l/len
  
  G1.p1 <- apply(data.training[1:G1.l,], 2, function(x) length(which(x == 1))/G1.l)
  G1.p2 <- apply(data.training[1:G1.l,], 2, function(x) length(which(x == -1))/G1.l)
  G1.p <- cbind(G1.p1, G1.p2)
  
  S.p1 <- apply(data.training[(G1.l+1):(G1.l+S.l),], 2, function(x) length(which(x == 1))/S.l)
  S.p2 <- apply(data.training[(G1.l+1):(G1.l+S.l),], 2, function(x) length(which(x == -1))/S.l)
  S.p <- cbind(S.p1, S.p2)
  
  G2M.p1 <- apply(data.training[(G1.l+S.l+1):len,], 2, function(x) length(which(x == 1))/G2M.l)
  G2M.p2 <- apply(data.training[(G1.l+S.l+1):len,], 2, function(x) length(which(x == -1))/G2M.l)
  G2M.p <- cbind(G2M.p1, G2M.p2)
  
  zero <- which(G1.p[,1] == 0 |G1.p[,1] == 1 | S.p[,1] == 0 |S.p[,1] == 1| G2M.p[,1] == 0| G2M.p[,1] == 1)
  if (length(zero) != 0){
    data.test <- data.test[,-zero]
    G1.p <- G1.p[-zero, ]
    S.p <- S.p[-zero, ]
    G2M.p <- G2M.p[-zero, ]
  }
  
  G1.score <- apply(data.test, 1, function(x) cellscore(x, G1.p, P.G1))
  S.score <- apply(data.test, 1, function(x) cellscore(x, S.p, P.S))
  G2M.score <- apply(data.test, 1, function(x) cellscore(x, G2M.p, P.G2M))
  
  result <- data.frame(G1.score = G1.score, S.score = S.score, G2M.score = G2M.score)
  
  return(result)
}