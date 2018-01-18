
# download.file("http://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz",
#               destfile = "sumstats-height.txt.gz")
# R.utils::gunzip("sumstats-height.txt.gz")

library(bigsnpr)

celiac <- snp_attach("backingfiles/celiacQC.rds")

sumstats <- data.table::fread("sumstats-height.txt", data.table = FALSE)

indSNP <- match(celiac$map$marker.ID, sumstats$MarkerName)
print(nrow(sumstats) - sum(!is.na(indSNP)))

same_ref

same <- same_ref(celiac$map$allele1, celiac$map$allele2, 
                 sumstats$Allele1[indSNP], sumstats$Allele2[indSNP])
sum(same, na.rm = TRUE)


G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
NCORES <- nb_cores()
betas <- sumstats$b[indSNP]
pval <- sumstats$p[indSNP]
lpS <- -log10(pval)


path <- "PRSice_linux.bugfix.20171116/"
# prefix <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel"
tmp <- tempfile()
snp_writeBed(celiac, paste0(tmp, ".bed"))

system.time(
  system(glue::glue("Rscript {path}PRSice.R",
                    " --prsice {path}PRSice_linux",
                    " --base sumstats-height.txt",
                    " --target {tmp}",
                    " --out {tmp}",
                    " --thread {NCORES}",
                    # " --full",
                    # " --fastscore",
                    # " --cov-file covar.txt",
                    # " --quantile 20",
                    " --pvalue p",
                    " --stat b",
                    " --beta T",
                    " --snp MarkerName",
                    " --A1 Allele1 --A2 Allele2",
                    " --binary-target T"))
)

file.prsice <- paste0(tmp, ".prsice")
cor.pval0 <- data.table::fread(file.prsice, select = "P")[["P"]]

library(dplyr)
thrs <- -log10(data.table::fread(file.prsice, select = "Threshold")[["Threshold"]])

system.time(
  ind.keep <- snp_clumping(G, infos.chr = CHR, S = lpS, thr.r2 = 0.1,
                           size = 250, is.size.in.bp = TRUE, infos.pos = POS,
                           exclude = which(is.na(same) | (lpS < min(thrs))), 
                           ncores = NCORES)
) # 2 -> 19, 1 -> 29

length(ind.keep)  # PRSice: 52946
system.time(
  prs <- snp_PRS(G, betas.keep = betas[ind.keep],
                 ind.test = rows_along(G),
                 same.keep = same[ind.keep], 
                 lpS.keep = lpS[ind.keep], 
                 ind.keep = ind.keep,
                 thr.list = thrs)
)

y01 <- celiac$fam$affection - 1
# cor.pval <- apply(prs, 2, function(scores) {
#   summary(glm(y01 ~ scores, family = "binomial"))$coefficients[2, 4]
# })
library(magrittr)
cor.pval2 <- prs %>%
  big_copy() %>%
  big_univLogReg(y01) %>%
  predict(log10 = FALSE)
# plot(cor.pval, cor.pval2)
# all.equal(cor.pval, cor.pval2)
all.equal(cor.pval2, cor.pval0)
plot(cor.pval2, cor.pval0, log = "xy"); abline(0, 1, col = "red")
cor(log10(cor.pval2), log10(cor.pval0))

CHR[lpS > 10]
plot(10^(-thrs), -log10(cor.pval2))


# obj.svd <- snp_autoSVD(G, CHR, POS, ncores = NCORES)
# write.table(cbind(celiac$fam[, 1:2], obj.svd$u), "covar.txt", 
#             quote = FALSE, row.names = FALSE)
# cor.pval2 <- prs %>%
#   big_copy() %>%
#   big_univLogReg(y01, covar.train = obj.svd$u) %>%
#   predict(log10 = FALSE)
# plot(10^(-thrs), -log10(cor.pval2))

# summary(glm(y01 ~ prs[, 1] + obj.svd$u, family = "binomial"))
