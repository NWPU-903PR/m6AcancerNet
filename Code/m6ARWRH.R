options(stringsAsFactors = FALSE)
library(Matrix)

genemat <- readRDS("./pcc_genemania.rds")
sitemat <- readRDS("./site_genemania.rds")
gene_site <- readRDS("./gene_site_genemania.rds")
M <- readRDS('./M_genemania.rds')

##sparse
M <- Matrix(data = M,sparse = TRUE)

#make p0
u0 <- Matrix(0,ncol = 1,nrow = dim(genemat)[1])
v0 <- Matrix(1,ncol = 1,nrow = dim(sitemat)[1])
v0 <- v0 / colSums(v0)
eta <- 0.5
geseeds <- c("METTL3","METTL14","WTAP","KIAA1429","RBM15","RBM15B","ZC3H13",
            "FTO","ALKBH5",
            "YTHDC1","YTHDC2","IGF2BP1","IGF2BP2","IGF2BP3","YTHDF1",
            "YTHDF2","YTHDF3","HNRNPA2B1","HNRNPC","RBMX")
geseeds <- geseeds[geseeds %in% rownames(genemat)]
topgenes <- matrix(0,ncol = length(geseeds),nrow = 100)

seed_connect <- function(gene_seed){
  loc <- match(gene_seed,rownames(genemat))
  seedrelation <- colnames(genemat)[genemat[loc,] == 1]
  return(seedrelation)}

#set regulator as source node
set_p0 <- function(geseed){
  u0[match(geseed,rownames(genemat)),] <- 1
  u0 <- u0 / colSums(u0)
  p0 <- rbind((1-eta)*u0,eta*v0)
  return(p0)
}

rwrh <- function(p0){
  gamma <- 0.5
  p1 <- p0
  p <- (1-gamma)*t(M)%*%p1+gamma*p0
  while(sum(abs(p-p1))>10^-10){
    p1 <- p
    p <-(1-gamma)*t(M)%*%p1+gamma*p0
  }
  p <- as.matrix(p)
  return(p)
}

top_gene <- function(p){
  u <- p[1:dim(u0)[1]]
  s <- sort(u,decreasing=TRUE,index.return=TRUE)
  top <- as.matrix(s$ix[1:100])
  topname <- rownames(genemat)[top]
  return(topname)
}


for(i in 1:length(geseeds)){
  p0 <- set_p0(geseeds[i])
  p <- rwrh(p0)
  topname <- top_gene(p)
  topgenes[,i] <- topname
}

#save result with txt
top_gene_100 <- unique(as.vector(topgenes))
con <- file('two_networks_topgenes.txt',open = "w")
for (i in 1:length(top_gene_100)){
  writeLines(top_gene_100[i],sep = '\n',con = con )
}
close(con)
