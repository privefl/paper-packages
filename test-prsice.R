
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
prefix <- "../thesis-celiac/Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel"

system.time(
  system(glue::glue("Rscript {path}PRSice.R",
                    " --prsice {path}PRSice_linux",
                    " --base sumstats-height.txt",
                    " --target {prefix}",
                    " --thread 6",
                    # " --full",
                    # " --fastscore",
                    " --pvalue p",
                    " --stat b",
                    " --beta T",
                    " --snp MarkerName",
                    " --A1 Allele1 --A2 Allele2",
                    " --binary-target T"))
)


cor.pval0 <- data.table::fread("PRSice.prsice", select = "P")[["P"]] 

library(dplyr)
thrs <- -log10(data.table::fread("PRSice.prsice", select = "Threshold")[["Threshold"]])

system.time(
  ind.keep <- snp_clumping(G, infos.chr = CHR, S = lpS, thr.r2 = 0.1,
                           size = 250, is.size.in.bp = TRUE, infos.pos = POS,
                           exclude = which(is.na(same) | (lpS < min(thrs))), 
                           ncores = NCORES)
)

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
cor.pval2 <- prs %>%
  big_copy() %>%
  big_univLogReg(y01) %>%
  predict(log10 = FALSE)
# plot(cor.pval, cor.pval2)
# all.equal(cor.pval, cor.pval2)
all.equal(cor.pval2, cor.pval0)
plot(cor.pval2, cor.pval0, log = "xy"); abline(0, 1, col = "red")


CHR[lpS > 10]
