library(bigsnpr)


# Write temporary PLINK files
tmp <- tempfile()
bed <- snp_writeBed(snp_attach("backingfiles/celiacQC.rds"), paste0(tmp, ".bed"))

library(SNPRelate)
# Transform files in our format and attach them
system.time(
  gds <- snpgdsBED2GDS(bed.fn = bed, 
                       fam.fn = sub("\\.bed$", ".fam", bed),
                       bim.fn = sub("\\.bed$", ".bim", bed),
                       out.gdsfn = paste0(tmp, ".gds"),
                       cvt.snpid = "int")
) # 33 sec

genofile <- snpgdsOpen(gds)
genofile
snpgdsClose(genofile)


# system.time(
#   rds <- snp_readBed(bed, tmp)
# ) # 43 sec

celiac <- snp_attach("backingfiles/celiacQC.rds")

G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
y01 <- celiac$fam$affection - 1
NCORES <- nb_cores()

# Pruning
system.time(
  snpset <- snpgdsLDpruning(genofile, method = "corr", ld.threshold = 0.2,
                            num.thread = NCORES)
) # 38 sec



ind.excl <- union(which(CHR == 2), snp_indLRLDR(CHR, POS))

system.time(
  ind.keep <- snp_pruning(G, CHR, size = 500, is.size.in.bp = TRUE,
                          exclude = snp_indLRLDR(CHR, POS),
                          thr.r2 = 0.1,
                          infos.pos = POS, ncores = NCORES)
) # 73 sec
length(ind.keep)

system.time(
  ind.keep2 <- snp_clumping(G, CHR, size = 500, is.size.in.bp = TRUE,
                            exclude = snp_indLRLDR(CHR, POS),
                            thr.r2 = 0.1,
                            infos.pos = POS, ncores = NCORES)
) # 66 sec
length(ind.keep2)

# system.time(
#   pca <- snpgdsPCA(genofile, snp.id = celiac$map$marker.ID[ind.keep],
#                    num.thread = NCORES, eigen.cnt = 10, algorithm = "randomized")
# ) # 622


system.time(
  obj.svd <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep, ncores = NCORES)
) # 148

system.time(
  obj.svd2 <- big_randomSVD(G, snp_scaleBinom(), ind.col = ind.keep2, ncores = NCORES)
) # 148

cowplot::plot_grid(
  plot(obj.svd), 
  plot(obj.svd2), 
  ncol = 1
)

# Get population from external files
pop.files <- list.files(path = "../Dubois2010_data/",
                        pattern = "cluster_*",
                        full.names = TRUE)
pop <- snp_getSampleInfos(celiac, pop.files)[[1]]
pop.names <- c("Netherlands", "Italy", "UK1", "UK2", "Finland")

library(ggplot2)
PCs <- 5:6
cowplot::plot_grid(
  plot(obj.svd, type = "scores", scores = PCs) +
    aes(color = pop.names[pop]) +
    labs(color = "Population", title = NULL),
  plot(obj.svd2, type = "scores", scores = PCs) +
    aes(color = pop.names[pop]) +
    labs(color = "Population", title = NULL), 
  ncol = 1
)

round(100 * cor(obj.svd$u, obj.svd2$u), 1)

round(100 * cor(pca$eigenvect, obj.svd$u), 1)

library(GWASTools)
(genofile <- GdsGenotypeReader(gds, YchromCode=24L, XYchromCode=25L))
scanID <- getScanID(genofile)
# sex <- getVariable(genofile, "sample.annot/sex")
phenotype <- getVariable(genofile, "sample.annot/phenotype")
covar <- obj.svd$u
names(covar) <- paste0("PC", 1:10)
scanAnnot <- ScanAnnotationDataFrame(
  cbind.data.frame(scanID, pheno = phenotype - 1)
)
snpID <- getSnpID(genofile)
chromosome <- getChromosome(genofile)
position <- getPosition(genofile)
alleleA <- getAlleleA(genofile)
alleleB <- getAlleleB(genofile)
snpAnnot <- SnpAnnotationDataFrame(data.frame(
  snpID = seq_along(snpID), chromosome, position, rsID = snpID,
  alleleA, alleleB, stringsAsFactors=FALSE),
  YchromCode=24L, XYchromCode=25L)

genoData <- GenotypeData(genofile, scanAnnot = scanAnnot, 
                         snpAnnot = snpAnnot)

system.time(
  gwas <- GWASTools::assocRegression(genoData, "pheno",
                                     model.type = "logistic",
                                     snpStart = 10e3)
) # 220 pour 5000 SNPs

gwas

system.time(
  obj.gwas <- big_univLogReg(G, y01, covar.train = obj.svd$u, 
                             ncores = NCORES)
) # 936.113


str(genoData)
