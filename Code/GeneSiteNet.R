library(org.Hs.eg.db)

# load m6Acomet Data
comet <- readRDS("dict_Bon_sm_trim.rds")
mgr <- readRDS("mgr.rds")

site_entrez <- unlist(mcols(mgr)$site_entrez,use.names = FALSE)
site_name <- mcols(mgr)$mod_name


entrezid <- keys(org.Hs.eg.db,keytype = "ENTREZID")
entrez_symbol <- AnnotationDbi::select(org.Hs.eg.db,keys =entrezid ,
                        columns = 'SYMBOL',keytype = 'ENTREZID')
site_symbol <- entrez_symbol[match(site_entrez,entrez_symbol$ENTREZID),]
rownames(site_symbol) <- site_name
site_symbol <- site_symbol[!is.na(site_symbol$ENTREZID),]

# load gene data
pcc <- readRDS('./pcc_genemania.rds')
gene_name <- rownames(pcc)

#site map to gene
site_map <- site_symbol[!is.na(match(site_symbol$SYMBOL,gene_name)),]
site_del <- site_symbol[is.na(match(site_symbol$SYMBOL,gene_name)),]

#check site
site_matrix <- readRDS('./site.rds')
del <- c()
check <- function(site){
  loc <- match(site,rownames(site_matrix))
  sum(site_matrix[loc,])
}
for (i in 1:dim(site_del)[1]){
  if (check(rownames(site_del)[i]) == 0){
    del <- c(del,rownames(site_del)[i])
  }
}
site_reserve <- site_symbol[-match(del,rownames(site_symbol)),]

# gene-site 
len_site <- length(rownames(site_reserve))
len_gene <- length(gene_name)
gene_site <- matrix(0,nrow = len_gene,ncol = len_site)
rownames(gene_site) <- gene_name
colnames(gene_site) <- rownames(site_reserve)

site <- rownames(site_map)
site_location <- match(site,rownames(site_reserve))
table(site_map$SYMBOL %in% gene_name)
gene_location <- match(site_map$SYMBOL,gene_name)
for (i in 1:dim(site_map)[1]){
  gene_site[gene_location[i],site_location[i]] <- 1
}

#check whether exists isolated genes
rowgene <- data.frame(name=rownames(pcc),sum=rowSums(pcc))
possible_del_gene <- rowgene[rowgene$sum == 0,]
del_gene <- possible_del_gene[!possible_del_gene$name %in% site_map$SYMBOL,]
gene_site <- gene_site[-match(del_gene$name,rownames(gene_site)),]

#delete  sites in site matrix
site_matrix <- site_matrix[-match(del,rownames(site_matrix)),-match(del,rownames(site_matrix))]

#delete genes in pcc
pcc <- pcc[-match(del_gene$name,rownames(pcc)),-match(del_gene$name,colnames(pcc))]

#save
saveRDS(pcc,'pcc_genemania.rds')
saveRDS(site_matrix,'site_genemania.rds')
saveRDS(gene_site,'gene_site_genemania.rds')
