#load data
Ag <- readRDS("./pcc_genemania.rds")
As <- readRDS("./site_genemania.rds")
B <- readRDS("./gene_site_genemania.rds")


colB <- colSums(B)
rowB <- rowSums(B)
lamda <- 0.5

# according row normalization
#Mg & Mgs
Mg <- matrix(nrow = dim(Ag)[1],ncol = dim(Ag)[2])
Mgs <- matrix(nrow = dim(B)[1],ncol = dim(B)[2])
rowAg <- rowSums(Ag)
for (i in 1:dim(Ag)[1]){
  if (rowB[i] == 0){
    Mg[i,] <- Ag[i,] / rowAg[i]
    Mgs[i,] <- 0
  }
  else{
    if (rowAg[i] != 0){
      Mg[i,] <- (1-lamda)*(Ag[i,] / rowAg[i])
      Mgs[i,] <- lamda * (B[i,] / rowB[i])
    }
    else{
      Mg[i,] <- 0
      Mgs[i,] <- B[i,] / rowB[i]
    }
  }
}

#Ms & Msg
Ms <- matrix(nrow = dim(As)[1],ncol = dim(As)[2])
Msg <- matrix(nrow = dim(B)[2],ncol = dim(B)[1])
rowAs <- rowSums(As)
for (i in 1:dim(As)[1]){
  if (colB[i] == 0){
    Ms[i,] <- As[i,] / rowAs[i]
    Msg[i,] <- 0
  }
  else{
    if(rowAs[i] == 0 ){
      Ms[i,] <- 0
      Msg[i,] <- B[,i] / colB[i]
    }
    else{
      Ms[i,] <- (1-lamda) * (As[i,] / rowAs[i])
      Msg[i,] <- lamda * (B[,i] / colB[i])
    }
  }
}

# transition matrix
M <- rbind(cbind(Mg,Mgs),cbind(Msg,Ms))

#check
rowM <- as.data.frame(rowSums(M))
colnames(rowM) <- 'sum'
rowM[is.nan(rowM$sum),]

saveRDS(M,'M_genemania.rds')
