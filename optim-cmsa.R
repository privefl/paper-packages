library(bigsnpr)

celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes
y01 <- celiac$fam$affection - 1
dim(G)

CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos

NCORES <- nb_cores()
ind.excl <- snp_indLRLDR(CHR, POS)
ind.keep <- snp_clumping(G, CHR, exclude = ind.excl, ncores = NCORES) 
obj.svd <- snp_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep, 
                         ncores = NCORES, verbose = TRUE) 

ind.train <- 1:12000
system.time(
  test <- big_spLogReg(G, y01[ind.train], ind.train = ind.train,
                       alpha = 0.999, dfmax = 20e3)
) # 490 sec

ind.test <- setdiff(rows_along(G), ind.train)
system.time(
  preds <- predict(test, G, ind.test)
) # 52 sec
system.time(
  aucs <- apply(preds, 2, function(pred) AUC(pred, y01[ind.test]))
)
summary(aucs[-1])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8156  0.8520  0.8596  0.8575  0.8725  0.8794 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.8156  0.8520  0.8596  0.8575  0.8725  0.8796 

plot(aucs[-1])

# b <- test$beta
# library(Matrix)
# rowSums(b != 0)
# colSums(b != 0)