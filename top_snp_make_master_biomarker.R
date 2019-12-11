#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

trait <- as.character(args[1])

library(data.table)
############ make master ######################
###############
#   chr1_22   #
###############

chromosomes <- c(1:22)
master1_22 <- c()
# make master files
for(ii in seq_along(chromosomes)){
  master.chr <- c()
  chr <- chromosomes[ii]
  z <- sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/z_files/top_snp_sex_combined_%s_chr%d_wim_1mb.z",
               trait, chr, trait, chr)
  ld <- sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.ld",
                trait, chr)
  bcor <- sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.bcor",
                  trait, chr)
  bdose <- sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.bdose",
                   trait, chr)
  bgen <- paste0("/fs/projects/ukbb/geno_v3/ukb_imp_chr", chr,"_v3.bgen")
  bgi <- paste0("/fs/projects/ukbb/geno_v3/ukb_imp_chr",chr,"_v3.bgen.bgi")
  sample <- rep("/fs/projects/ukbb/geno_v3/ukb22627_imp_chr1-22_v3_s487395.sample",length(chr))
  incl <- rep(sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/sex_combined_BOLT_%s.incl",trait), length(chr))
  total_records <- as.integer(system2("wc",
                                      args = c("-l",
                                               incl,
                                               " | awk '{print $1}'"),
                                      stdout = TRUE))
  n_samples <- rep(total_records, length(chr))
  if(file.exists(z)){
    n.z <- nrow(read.table(z, header = TRUE))
    if(n.z > 1){
      master.chr <- cbind(z,ld,bcor,bgen,bgi,bdose,sample,incl,n_samples)
      master1_22 <- rbind(master1_22, master.chr)
    }
  }
}

###########
#   chrX  #
###########

##########  non-PAR #############

chromosomes <- 23

z <- sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/z_files/top_snp_sex_combined_%s_chr%d_wim_1mb.z",
             trait, chromosomes,trait, chromosomes)
ld <- sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.ld",
              trait, chromosomes)
bcor <- sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.bcor",
                trait, chromosomes)
bdose <- sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.bdose",
                 trait, chromosomes)
bgen <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrX_v3.bgen"
bgi <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrX_v3.bgen.bgi"
sample <- "/fs/projects/ukbb/geno_v3/ukb22627_imp_chrX_v3_s486743.sample"
incl <- rep(sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/sex_combined_BOLT_%s.incl",trait), length(chromosomes))
n_samples <- rep(total_records, length(chromosomes))

if(file.exists(z)){
  n.z <- nrow(read.table(z, header = TRUE))
  if(n.z > 1){
    masterX <- cbind(z,ld,bcor,bgen,bgi,bdose,sample,incl,n_samples)
    master <- rbind(master1_22, masterX)
  }
} else {
  master <- master1_22
}

write.table(master,file=paste0("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/",trait,"/sex_combined/chr1-23_master"), quote = FALSE, sep = ";", row.names = FALSE, col.names = TRUE)


########## PAR ##########
#chromosomes <- "PAR"
#z <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/z_files/top_snp_sex_combined_%s_%s_wim_1mb.z",
#             trait, chromosomes,trait, chromosomes)
#
#ld <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/ld/top_snp_sex_combined_win_1mb.ld",
#              trait, chromosomes)
#bcor <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/ld/top_snp_sex_combined_win_1mb.bcor",
#                trait, chromosomes)
#bdose <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/%s/ld/top_snp_sex_combined_win_1mb.bdose",
#                 trait, chromosomes)
#bgen <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrXY_v3.bgen"
#bgi <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrXY_v3.bgen.bgi"
#sample <- "/fs/projects/ukbb/geno_v3/ukb22627_imp_chrXY_v3_s486429.sample"
#incl <- rep(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/sex_combined_BOLT_%s.incl",trait), length(chromosomes))
#n_samples <- rep(total_records, length(chromosomes))
#
#masterPAR <- cbind(z,ld,bcor,bgen,bgi,bdose,sample,incl,n_samples)
#


