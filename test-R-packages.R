library(bigsnpr)

# Write temporary PLINK files
tmp <- tempfile()
bed <- snp_writeBed(snp_attach("backingfiles/celiacQC.rds"), paste0(tmp, ".bed"))

# gaston: Genetic Data Handling (QC, GRM, LD, PCA) & Linear Mixed Models
library(gaston) 

system.time(
  x <- read.bed.matrix(bed)
)

x[1, 1]


library(snpStats)
snp <- read.plink(bed, 
                  sub("\\.bed$", ".fam", bed),
                  sub("\\.bed$", ".bim", bed)) ## All in memory


library(GenABEL)  # pas updated depuis 2013

?GenABEL::export.plink

plink <- download_plink()
system(glue::glue("{plink} --bfile {tmp} --recode --out {tmp}"))
file.size(paste0(tmp, ".ped")) / 1024^3
system.time(
  test <- convert.snp.ped(pedfile = paste0(tmp, ".ped"),
                          mapfile = paste0(tmp, ".map"),
                          out = "test.dat"))
)
