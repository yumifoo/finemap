# make input file master file for ld store of each region
# i.e. the master file will be a two lines file, one header, one content


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
trait <- as.character(args[1])
maf <- 0.01


library(data.table)
source("/fs/projects/ukbb/yu/src/BOLT2z.R")

#########################
###     autosomes     ###
#########################

chromosomes <- 1:22

##### load GWAS results ######
res <- fread(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/BOLT_sex_combined_s343695_WhiteBritish.%s.bgen.stats.gz", trait))
res$MAF <- ifelse(res$A1FREQ < 0.5, res$A1FREQ, 1-res$A1FREQ)
res <- res[res$MAF > maf, ]

####### make z files for each region ############
for (ii in seq_along(chromosomes)){
  chr <- chromosomes[ii]
  region.file <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/regions/merged_region_1mb.txt",trait, chr)
  if(file.exists(region.file)){
    n.region <- as.integer(system2("wc",
                                   args = c("-l",
                                            region.file,
                                            " | awk '{print $1}'"),
                                   stdout = TRUE))
    if(n.region > 0){
      regions <- read.table(region.file, header = FALSE)
      auto.res <- subset(res, CHR == chr)
      filenames <- c()
      infile.ld <- c()
      infile.finemap <- c()
      
      for (jj in 1:nrow(regions)){
        subres <- subset(auto.res, BP > regions[jj,1] & BP < regions[jj,2])
        output <- sprintf("BOLT_sex_combined_s343695_WhiteBritish_MAF_0.01_%s_chr%d_%d_%d",
                          trait, chr,regions[jj,1],regions[jj,2])
        BOLT2z(subres, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/z_files/%s",
                               trait,chr,output))
        filenames <- c(filenames,output)
        print(paste("The number of SNPs in", chr, regions[jj,1], regions[jj,2], ":", nrow(subres)))
        
        
        ####### make master files for each region ########
        
        z <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/z_files/%s.z",
                     trait,chr,output)
        #ld <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/%s.ld",
        #              trait,chr,output)
        bcor <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/%s.bcor",
                        trait,chr,output)
        bdose <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/%s.bdose",
                         trait,chr,output)
        bgen <- sprintf("/fs/projects/ukbb/geno_v3/ukb_imp_chr%d_v3.bgen",chr)
        bgi <- sprintf("/fs/projects/ukbb/geno_v3/ukb_imp_chr%d_v3.bgen.bgi",chr)
        sample <- "/fs/projects/ukbb/geno_v3/ukb22627_imp_chr1-22_v3_s487395.sample"
        incl <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/sex_combined_BOLT_%s.incl",trait)
        n_samples <- as.integer(system2("wc",
                                        args = c("-l",
                                                 incl,
                                                 " | awk '{print $1}'"),
                                        stdout = TRUE))
        snp <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/%s.snp",
                       trait,chr,output)
        config <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/%s.config",
                          trait,chr,output)
        cred <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/%s.cred",
                        trait,chr,output)
        log <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/%s.log",
                       trait,chr,output)
        master <- cbind(z,bcor, bdose,bgen,bgi,sample,incl,n_samples,snp,config,cred,log)
        write.table(master, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master/master_chr%d_%d_%d",
                                           trait, chr, chr, regions[jj,1], regions[jj,2]),
                    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ";")
        #infile[jj] <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master/master_chr%d_%d_%d",
        #                      trait, chr, chr, regions[jj,1], regions[jj,2])
        #write.table(infile, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master_list.txt",trait,chr),
        #            quote = FALSE, row.names = FALSE, col.names = FALSE)
        ########### make ld commands ##############
        command_ld <- sprintf("/fs/projects/ukbb/christian/binaries/150419/ldstore_v2.0b_x86_64 --in-files /fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master/master_chr%d_%d_%d --write-bcor --write-bdose --n-threads 50 --cpu-mem 100",
                              trait, chr, chr,regions[jj,1], regions[jj,2])
        write.table(command_ld,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/ldstore_chr%d_%d_%d.sh",
                                       trait, chr, chr, regions[jj,1], regions[jj,2]), 
                    quote = FALSE, row.names = FALSE, col.names = "#!/bin/bash")
        infile.ld[jj] <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/ldstore_chr%d_%d_%d.sh",
                                 trait, chr, chr, regions[jj,1], regions[jj,2])
        write.table(infile.ld, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/ldstore_list.txt",trait,chr),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        ########### make finemap commands ##############
        command_finemap <- sprintf("/fs/projects/ukbb/christian/binaries/180419/finemap_v1.4_x86_64 --cond --in-files /fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master/master_chr%d_%d_%d --log --n-causal-snps 30", 
                                   trait, chr, chr,regions[jj,1], regions[jj,2])
        write.table(command_finemap, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/finemap_cond_chr%d_%d_%d.sh",
                                             trait, chr, chr, regions[jj,1], regions[jj,2]),
                    quote = FALSE, row.names = FALSE, col.names = "#!/bin/bash")
        infile.finemap[jj] <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/finemap_cond_chr%d_%d_%d.sh",
                                      trait, chr, chr, regions[jj,1], regions[jj,2])
        write.table(infile.finemap, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/finemap_cond_list.txt",trait,chr),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
    }
  }
}
  
  
  ###########################
  ###     chromsome X     ###
  ###########################
  
  chr <- 23
  non.PAR.list <- fread("/fs/projects/ukbb/yu/ukb_imp_chrX_snplist.txt")
  region.file <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/regions/merged_region_1mb.txt",trait, chr)
  if(file.exists(region.file)){
    n.region <- as.integer(system2("wc",
                                   args = c("-l",
                                            region.file,
                                            " | awk '{print $1}'"),
                                   stdout = TRUE))
    if(n.region > 0){
      regions <- read.table(region.file, header = FALSE)
      res.x <- subset(res, CHR == chr)
      auto.res <- subset(res.x, BP %in% non.PAR.list$position)
      filenames <- c()
      infile.ld <- c()  
      infile.finemap <- c()  
      
      for (jj in 1:nrow(regions)){
        subres <- subset(auto.res, BP > regions[jj,1] & BP < regions[jj,2])
        output <- sprintf("BOLT_sex_combined_s343695_WhiteBritish_MAF_0.01_%s_chr%d_%d_%d",
                          trait, chr,regions[jj,1],regions[jj,2])
        BOLT2z(subres, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/z_files/%s",
                               trait,chr,output))
        filenames <- c(filenames,output)
        print(paste("The number of SNPs in", chr, regions[jj,1], regions[jj,2], ":", nrow(subres)))
        
        
        ####### make master files for each region ########
        
        z <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/z_files/%s.z",
                     trait,chr,output)
        #ld <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/%s.ld",
        #trait,chr,output)
        bcor <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/%s.bcor",
                        trait,chr,output)
        bdose <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/%s.bdose",
                         trait,chr,output)
        bgen <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrX_v3.bgen"
        bgi <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrX_v3.bgen.bgi"
        sample <- "/fs/projects/ukbb/geno_v3/ukb22627_imp_chrX_v3_s486743.sample"
        incl <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/sex_combined_BOLT_%s.incl",trait)
        n_samples <- as.integer(system2("wc",
                                        args = c("-l",
                                                 incl,
                                                 " | awk '{print $1}'"),
                                        stdout = TRUE))
        snp <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/%s.snp",
                       trait,chr,output)
        config <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/%s.config",
                          trait,chr,output)
        cred <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/%s.cred",
                        trait,chr,output)
        log <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/%s.log",
                       trait,chr,output)
        master <- cbind(z,bcor, bdose,bgen,bgi,sample,incl,n_samples,snp,config,cred,log)
        write.table(master, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master/master_chr%d_%d_%d",
                                           trait, chr, chr, regions[jj,1], regions[jj,2]),
                    quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ";")
        #infile[jj] <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master/master_chr%d_%d_%d",
        #                      trait, chr, chr, regions[jj,1], regions[jj,2])
        #write.table(infile, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master_list.txt",trait,chr),
        #            quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        ########### make ld commands ##############
        command_ld <- sprintf("/fs/projects/ukbb/christian/binaries/150419/ldstore_v2.0b_x86_64 --in-files /fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master/master_chr%d_%d_%d --write-bcor --write-bdose --n-threads 50 --cpu-mem 100",
                              trait, chr, chr,regions[jj,1], regions[jj,2])
        write.table(command_ld,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/ldstore_chr%d_%d_%d.sh",
                                       trait, chr, chr, regions[jj,1], regions[jj,2]), 
                    quote = FALSE, row.names = FALSE, col.names = "#!/bin/bash")
        infile.ld[jj] <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/ldstore_chr%d_%d_%d.sh",
                                 trait, chr, chr, regions[jj,1], regions[jj,2])
        write.table(infile.ld, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/ldstore_list.txt",trait,chr),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
        
        ########### make finemap commands ##############
        command_finemap <- sprintf("/fs/projects/ukbb/christian/binaries/180419/finemap_v1.4_x86_64 --cond --in-files /fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master/master_chr%d_%d_%d --log --n-causal-snps 30", 
                                   trait, chr, chr,regions[jj,1], regions[jj,2])
        write.table(command_finemap, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/finemap_cond_chr%d_%d_%d.sh",
                                             trait, chr, chr, regions[jj,1], regions[jj,2]),
                    quote = FALSE, row.names = FALSE, col.names = "#!/bin/bash")
        infile.finemap[jj] <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/finemap_cond_chr%d_%d_%d.sh",
                                      trait, chr, chr, regions[jj,1], regions[jj,2])
        write.table(infile.finemap, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/finemap_cond_list.txt",trait,chr),
                    quote = FALSE, row.names = FALSE, col.names = FALSE)
      }
    }
}


###########################
###        PAR          ###
###########################

#chr <- "PAR"
#PAR.list <- fread("/fs/projects/ukbb/yu/ukb_imp_chrXY_snplist.txt")
#
#regions <- read.table(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/regions/merged_region_1mb.txt",
#                              trait, chr), header = FALSE)
#auto.res <- subset(res.x, BP %in% PAR.list$position)
#filenames <- c()
#infile.ld <- c()  
#infile.finemap <- c()  
#
#for (jj in 1:nrow(regions)){
#  subres <- subset(auto.res, BP > regions[jj,1] & BP < regions[jj,2])
#  output <- sprintf("BOLT_sex_combined_s343695_WhiteBritish_MAF_0.01_%s_%s_%d_%d",
#                    trait, chr,regions[jj,1],regions[jj,2])
#  BOLT2z(subres, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/z_files/%s",
#                         trait,chr,output))
#  filenames <- c(filenames,output)
#  print(paste("The number of SNPs in", chr, regions[jj,1], regions[jj,2], ":", nrow(subres)))
  
  
  ####### make master files for each region ########
  
#  z <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/z_files/%s.z",
#               trait,chr,output)
  #ld <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/ld/%s.ld",
  #trait,chr,output)
#  bcor <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/ld/%s.bcor",
#                  trait,chr,output)
#  bdose <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/ld/%s.bdose",
#                   trait,chr,output)
#  bgen <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrXY_v3.bgen"
#  bgi <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrXY_v3.bgen.bgi"
#  sample <- "/fs/projects/ukbb/geno_v3/ukb22627_imp_chrXY_v3_s486429.sample"
#  incl <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/sex_combined_BOLT_%s.incl",trait)
#  n_samples <- as.integer(system2("wc",
#                                  args = c("-l",
#                                           incl,
#                                           " | awk '{print $1}'"),
#                                  stdout = TRUE))
#  snp <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/finemap/%s.snp",
#                 trait,chr,output)
#  config <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/finemap/%s.config",
#                    trait,chr,output)
#  cred <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/finemap/%s.cred",
#                  trait,chr,output)
#  log <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/finemap/%s.log",
#                 trait,chr,output)
#  master <- cbind(z,bcor, bdose,bgen,bgi,sample,incl,n_samples,snp,config,cred,log)
#  write.table(master, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/master/master_%s_%d_%d",
#                                     trait, chr, chr, regions[jj,1], regions[jj,2]),
#              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ";")
  #infile[jj] <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master/master_chr%d_%d_%d",
  #                      trait, chr, chr, regions[jj,1], regions[jj,2])
  #write.table(infile, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/master_list.txt",trait,chr),
  #            quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ########### make ld commands ##############
#  command_ld <- sprintf("/fs/projects/ukbb/christian/binaries/150419/ldstore_v2.0b_x86_64 --in-files /fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/master/master_%s_%d_%d --write-bcor --write-bdose --n-threads 50 --cpu-mem 100",
#                        trait, chr, chr,regions[jj,1], regions[jj,2])
#  write.table(command_ld,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/ld/ldstore_%s_%d_%d.sh",
#                                 trait, chr, chr, regions[jj,1], regions[jj,2]), 
#              quote = FALSE, row.names = FALSE, col.names = "#!/bin/bash")
#  infile.ld[jj] <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/ld/ldstore_%s_%d_%d.sh",
#                           trait, chr, chr, regions[jj,1], regions[jj,2])
#  write.table(infile.ld, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/ld/ldstore_list.txt",trait,chr),
#              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  ########### make finemap commands ##############
#  command_finemap <- sprintf("/fs/projects/ukbb/christian/binaries/180419/finemap_v1.4_x86_64 --cond --in-files /fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/master/master_%s_%d_%d --log --n-causal-snps 30", 
#                             trait, chr, chr,regions[jj,1], regions[jj,2])
#  write.table(command_finemap, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/finemap/finemap_cond_%s_%d_%d.sh",
#                                       trait, chr, chr, regions[jj,1], regions[jj,2]),
#              quote = FALSE, row.names = FALSE, col.names = "#!/bin/bash")
#  infile.finemap[jj] <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/finemap/finemap_cond_%s_%d_%d.sh",
#                                trait, chr, chr, regions[jj,1], regions[jj,2])
#  write.table(infile.finemap, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/finemap/finemap_cond_list.txt",trait,chr),
#              quote = FALSE, row.names = FALSE, col.names = FALSE)
#}





