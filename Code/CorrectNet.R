## Loding reference gene functional interaction network data
pcc.bin <- readRDS('./pcc_bin.rds')
load('./GeneMania.rda')
genename.pcc <- rownames(pcc.bin)
genename.ppi <- rownames(originalppi)
cover <- genename.pcc[genename.pcc %in% genename.ppi]
pcc.bin <- pcc.bin[match(cover, genename.pcc), match(cover, genename.pcc)]
originalppi <- originalppi[match(cover, genename.ppi), match(cover, genename.ppi)]
pcc.adj <- pcc.bin * originalppi
