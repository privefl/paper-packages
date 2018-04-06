library(glue)
library(foreach)

# http://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html
hapgen2 <- "~/Téléchargements/hapgen2_x86_64/hapgen2"
# https://mathgen.stats.ox.ac.uk/wtccc-software/HM3.tgz
path <- file.path(dirname(hapgen2), "CEU.0908.impute.files")
# http://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html
gtool <- "~/Téléchargements/hapgen2_x86_64/gtool_v0.7.5_x86_64/gtool"
# https://www.cog-genomics.org/plink2
plink <- bigsnpr::download_plink()
simus.path <- "simus_hapgen2"
dir.create(simus.path)

get_simu <- function(i, n = 1000) {
  system <- function(command) base::system(command, ignore.stdout = TRUE)
  # prefix of created files
  prefix <- glue("{simus.path}/simu{i}.out")
  # simulate with hapgen2 in Oxford format
  cmd.hapgen2 <- glue(
    "{hapgen2}",
    " -m {path}/genetic_map_chr2_combined_b36.txt",
    " -l {path}/CEU.0908.chr2.legend",
    " -h {path}/CEU.0908.chr2.hap",
    " -o {prefix}.gz",
    " -dl 205 1 1 1",
    " -n {n} 1", 
    " -no_haps_output"
  )
  system(cmd.hapgen2)
  # convert to ped.gz/map.gz
  cmd.gtool <- glue(
    "{gtool} -G",
    " --g {prefix}.controls.gen.gz", 
    " --s {prefix}.controls.sample"
  )
  system(cmd.gtool)
  # uncompress
  system(glue("gzip -d {prefix}.controls.gen.map.gz"))
  system(glue("gzip -d {prefix}.controls.gen.ped.gz"))
  # convert to bed/bim/fam
  system(glue("{plink} --file {prefix}.controls.gen --make-bed --out {prefix}"))
  # cleaning all but bed/bim/fam
  # https://stackoverflow.com/a/23179318/6103040
  list <- grep(pattern = "^[^.]+$|\\.(?!(bed|bim|fam)$)([^.]+$)", 
               x = list.files(simus.path, full.names = TRUE), 
               perl = TRUE, value = TRUE)
  file.remove(grep(basename(prefix), list, value = TRUE, fixed = TRUE))
}

get_simu(2, 200)

doParallel::registerDoParallel(parallel::detectCores() / 2)
foreach(ic = 67) %dopar% {
  # do not launch simultaneous jobs, cause of seeds on second
  while (round(lubridate::second(Sys.time())) != ((3 * ic) %% 60)) {
    Sys.sleep(0.5)
  }
  # launch job
  get_simu(ic)
  NULL
}
doParallel::stopImplicitCluster()

### Merge

# change names of individuals
for (i in 1:500) {
  print(i)
  famfile <- glue("{simus.path}/simu{i}.out.fam")
  df <- data.table::fread(famfile, data.table = FALSE)
  df[[1]] <- glue("{df[[1]]}_{i}")
  df[[2]] <- glue("{df[[2]]}_{i}")
  bigsnpr:::write.table2(df, famfile)
}


library(bigsnpr)

merge_plink <- function(n) {
  write(glue("{simus.path}/simu{2:n}.out"), tmpfile <- tempfile())
  system(glue(
    "{plink}",
    " --bfile {simus.path}/simu1.out",
    " --merge-list {tmpfile}",
    " --make-bed", 
    " --out {simus.path}/simus_1_{n}",
    " --maf 0.01"
  ))
}
merge_plink(10)
merge_plink(30)
merge_plink(100)
merge_plink(500)
prefix <- glue("{simus.path}/simu")

file.remove(paste0("backingfiles/test", c(".bk", ".rds")))
snp_readBed("simus_hapgen2/simus_1_10.bed", "test")
test <- snp_attach("backingfiles/test.rds")
big_counts(test$genotypes, ind.col = 1:10)
