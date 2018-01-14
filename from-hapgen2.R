library(bigsnpr)

rds <- snp_readBed("simus_hapgen2/simus_1_10.bed",
                   tmp <- tempfile())

bigsnp <- snp_attach(rds)
G <- bigsnp$genotypes
dim(G)
CHR <- bigsnp$map$chromosome
rle(CHR)
POS <- bigsnp$map$physical.pos

counts <- big_counts(G)
counts[, 1:10]
sum(counts[4, ])
which(counts[4, ] != 0)
indNA <- which(counts[4, ] != 0)

NCORES <- nb_cores()
ind.rm <- snp_indLRLDR(CHR, POS)
ind.keep <- snp_clumping(G, CHR, exclude = union(indNA, ind.rm), ncores = NCORES)
