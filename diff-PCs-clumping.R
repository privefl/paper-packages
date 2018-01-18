library(bigsnpr)
bigsnp <- snp_attach("backingfiles/dogs.rds")
G <- bigsnp$genotypes
CHR <- bigsnp$map$chromosome
POS <- bigsnp$map$physical.pos
NCORES <- nb_cores()

system.time(
  ind.keep <- snp_pruning(G, CHR, size = 500, is.size.in.bp = TRUE,
                          # thr.r2 = 0.1,
                          infos.pos = POS, ncores = NCORES)
) # 9 sec
length(ind.keep)

system.time(
  ind.keep2 <- snp_clumping(G, CHR, size = 500, is.size.in.bp = TRUE,
                            # thr.r2 = 0.1,
                            infos.pos = POS, ncores = NCORES)
) # 6 sec
length(ind.keep2)


system.time(
  obj.svd <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep, k = 50, ncores = NCORES)
) # 8 / 18

system.time(
  obj.svd2 <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep2, k = 50, ncores = NCORES)
) # 7 / 18

diag(round(100 * cor(obj.svd$u, obj.svd2$u), 1))

cowplot::plot_grid(
  plot(obj.svd), 
  plot(obj.svd2), 
  ncol = 1
)


library(ggplot2)
PCs <- 26:27
cowplot::plot_grid(
  plot(obj.svd, type = "scores", scores = PCs),
  plot(obj.svd2, type = "scores", scores = PCs), 
  ncol = 1
)

