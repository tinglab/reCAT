# in this function, you should input a matrix that you want to calculate bayes_score and mean_score.

# the matrix is transposition of test_exp(you use in get the ordIndex)

# ola_mES_2i.RData in R is a sample

  
source("bayes_score_function.R")
source("mean_score_function.R")
load("../data/tranning.RData")
load("../data/ccgene.RData")

get_score <- function(talentDE)
{
    rname = row.names(talentDE)
    if ((length(grep("ensg", tolower(rname[1]))) == 0) & (length(grep("ensmusg", tolower(rname[1]))) == 0))
    {
        rname = geneReferTable[match(tolower(rname), tolower(geneReferTable[,4])),3]
    }else
    {
        if (length(grep("ensg", tolower(rname[1]))) != 0)
        {
            rname = geneReferTable[match(tolower(rname), tolower(geneReferTable[,1])),3]
        }
    }
    row.names(talentDE) = rname

    y <- rep(0,length(training.data[1,]))
    cname <- colnames(training.data)
    G1.id <- grep("G1", cname)
    S.id <- grep("S", cname)
    G2M.id <- grep("G2", cname)
    y[G1.id] <- 1
    y[S.id] <- 2
    y[G2M.id] <- 3
    
    markers <- getPairs(talentDE, G1.marker.pairs, 
                        S.marker.pairs, G2M.marker.pairs, genes.training)
    vector.training <- t(getVector(training.data, markers))
    
    vector.data <- t(getVector(talentDE, markers))
    
    pred <- Naive(vector.training, length(G1.id), length(S.id), length(G2M.id), dim(training.data)[2], vector.data)
    
    result <- averageScoreall(talentDE, allgene)
    
    return(list(bayes_score = pred, mean_score = result))
}