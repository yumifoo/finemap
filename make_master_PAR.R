#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
trait <- as.character(args[1])
maf <- 0.01

library(data.table)

regions <- read.table(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/regions/merged_region_1mb.txt",
                              trait), header = FALSE)
##### load GWAS results ######
res <- fread(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/BOLT_sex_combined_s343695_WhiteBritish.%s.bgen.stats.gz", trait))
res$MAF <- ifelse(res$A1FREQ < 0.5, res$A1FREQ, 1-res$A1FREQ)
res <- res[res$MAF > maf, ]

par.list <- fread("/fs/projects/ukbb/yu/ukb_imp_chrXY_snplist.txt")

res.par <- subset(res, SNP %in% par.list$rsid)
res.par$CHR <- rep("XY",nrow(res.par))

for(jj in 1:nrow(regions)){
  subres <- subset(res.par, BP > regions[jj,1] & BP < regions[jj,2])
  output <- sprintf("BOLT_sex_combined_s343695_WhiteBritish_MAF_0.01_%s_PAR_%d_%d",
                    trait, regions[jj,1],regions[jj,2])
  z <- data.frame(rsid = subres$SNP,
                  chromosome = subres$CHR,
                  position = subres$BP,
                  allele1 = subres$ALLELE1,
                  allele2 = subres$ALLELE0,
                  maf = subres$MAF,
                  beta = subres$BETA,
                  se = subres$SE)
  write.table(z, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/z_files/%s.z",
                         trait,output), 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  ####### make master files for each region ########
  
  z <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/z_files/%s.z",
               trait,output)
  #ld <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/%s.ld",
  #trait,chr,output)
  bcor <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/ld/%s.bcor",
                  trait,output)
  bdose <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/ld/%s.bdose",
                   trait,output)
  bgen <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrXY_v3.bgen"
  bgi <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrXY_v3.bgen.bgi"
  sample <- "/fs/projects/ukbb/geno_v3/ukb22627_imp_chrXY_v3_s486429.sample"
  incl <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/sex_combined_BOLT_PAR_%s.incl",trait)
  n_samples <- as.integer(system2("wc",
                                  args = c("-l",
                                           incl,
                                           " | awk '{print $1}'"),
                                  stdout = TRUE))
  snp <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/finemap/%s.snp",
                 trait,output)
  config <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/finemap/%s.config",
                    trait,output)
  cred <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/finemap/%s.cred",
                  trait,output)
  log <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/finemap/%s.log",
                 trait,output)
  master <- cbind(z,bcor, bdose,bgen,bgi,sample,incl,n_samples,snp,config,cred,log)
  write.table(master, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/master/master_PAR_%d_%d",
                                     trait, regions[jj,1], regions[jj,2]),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ";")
  ############ make ld commands ################
  command_ld <- sprintf("/fs/projects/ukbb/christian/binaries/150419/ldstore_v2.0b_x86_64 --in-files /fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/master/master_PAR_%d_%d --write-bcor --write-bdose --n-threads 50 --cpu-mem 100",
                        trait, regions[jj,1], regions[jj,2])
  write.table(command_ld,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/ld/ldstore_PAR_%d_%d.sh",
                                 trait, regions[jj,1], regions[jj,2]), 
              quote = FALSE, row.names = FALSE, col.names = "#!/bin/bash")
  ########### make finemap commands #############
  command_finemap <- sprintf("/fs/projects/ukbb/christian/binaries/180419/finemap_v1.4_x86_64 --cond --in-files /fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/master/master_PAR_%d_%d --log --n-causal-snps 30", 
                             trait, regions[jj,1], regions[jj,2])
  write.table(command_finemap, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/finemap/finemap_cond_PAR_%d_%d.sh",
                                       trait, regions[jj,1], regions[jj,2]),
              quote = FALSE, row.names = FALSE, col.names = "#!/bin/bash")
                             
}

