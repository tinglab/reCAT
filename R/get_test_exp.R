ensIDToName_human <- as.matrix(read.table("../ensIDToName(Human).txt",header=FALSE,sep="\t",comment="",quote=""))
ensIDToName_mouse <- as.matrix(read.table("../ensIDToName(Mouse).txt",header=FALSE,sep="\t",comment="",quote=""))
hcyclebaseGeneList <- as.matrix(read.table("../cyclebaseGeneList",header=T,sep="\t",comment="",quote=""))[,1]
cbGeneIdx <- (match(tolower(ensIDToName_human[,2]),tolower(ensIDToName_mouse[,2])))
geneReferTable <- cbind(ensIDToName_human,ensIDToName_mouse[cbGeneIdx,])
hcyclebaseGeneList2 <- geneReferTable[match(tolower(hcyclebaseGeneList), tolower(geneReferTable[,4])),1]
hcyclebaseGeneList3 <- geneReferTable[match(tolower(hcyclebaseGeneList), tolower(geneReferTable[,4])),3]

get_test_exp <- function(data)
{
    rname = tolower(row.names(data))
    if ((length(grep("ensg", rname[1])) == 0) & (length(grep("ensmusg", rname[1])) == 0))
    {
        id1 = match(tolower(hcyclebaseGeneList), rname)
        id1 = id1[which(!is.na(id1))]
        rname = tolower(row.names(data[id1,]))
        id2 = match(rname, tolower(geneReferTable[,4]))
        id2 = which(!is.na(id2))
        id = id1[id2]
        return(t(data[id,]))
    }else
    {
        if (length(grep("ensg", rname[1])) != 0)
        {
            id1 = match(tolower(hcyclebaseGeneList2), rname)
            id1 = id1[which(!is.na(id1))]
            rname = tolower(row.names(data[id1,]))
            id2 = match(rname, tolower(geneReferTable[,1]))
            id2 = which(!is.na(id2))
            id = id1[id2]
            return(t(data[id,]))
        }else
        {
            id1 = match(tolower(hcyclebaseGeneList3), rname)
            id1 = id1[which(!is.na(id1))]
            rname = tolower(row.names(data[id1,]))
            id2 = match(rname, tolower(geneReferTable[,3]))
            id2 = which(!is.na(id2))
            id = id1[id2]
            return(t(data[id,]))
        }
    }
}