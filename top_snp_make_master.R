#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

trait <- as.character(args[1])

library(data.table)
############ make master ######################
###############
#   chr1_22   #
###############

chromosomes <- c(1:22)

# make master files

z <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/z_files/top_snp_sex_combined_%s_chr%d_wim_1mb.z",
             trait, chromosomes, trait, chromosomes)
ld <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.ld",
              trait, chromosomes)
bcor <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.bcor",
                trait, chromosomes)
bdose <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.bdose",
                 trait, chromosomes)
bgen <- paste0("/fs/projects/ukbb/geno_v3/ukb_imp_chr", chromosomes,"_v3.bgen")
bgi <- paste0("/fs/projects/ukbb/geno_v3/ukb_imp_chr",chromosomes,"_v3.bgen.bgi")
sample <- rep("/fs/projects/ukbb/geno_v3/ukb22627_imp_chr1-22_v3_s487395.sample",length(chromosomes))
incl <- rep(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/sex_combined_BOLT_%s.incl",trait), length(chromosomes))
target_file   <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/sex_combined_BOLT_%s.incl",trait)
total_records <- as.integer(system2("wc",
                                    args = c("-l",
                                             target_file,
                                             " | awk '{print $1}'"),
                                    stdout = TRUE))

n_samples <- rep(total_records, length(chromosomes))

master1_22 <- cbind(z,ld,bcor,bgen,bgi,bdose,sample,incl,n_samples)

###########
#   chrX  #
###########

##########  non-PAR #############

chromosomes <- 23

z <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/z_files/top_snp_sex_combined_%s_chr%d_wim_1mb.z",
             trait, chromosomes,trait, chromosomes)
ld <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.ld",
              trait, chromosomes)
bcor <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.bcor",
                trait, chromosomes)
bdose <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.bdose",
                 trait, chromosomes)
bgen <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrX_v3.bgen"
bgi <- "/fs/projects/ukbb/geno_v3/ukb_imp_chrX_v3.bgen.bgi"
sample <- "/fs/projects/ukbb/geno_v3/ukb22627_imp_chrX_v3_s486743.sample"
incl <- rep(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/sex_combined_BOLT_%s.incl",trait), length(chromosomes))
n_samples <- rep(total_records, length(chromosomes))

masterX <- cbind(z,ld,bcor,bgen,bgi,bdose,sample,incl,n_samples)

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
master <- rbind(master1_22, masterX)
#
write.table(master,file=paste0("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/",trait,"/sex_combined/chr1-23_master"), quote = FALSE, sep = ";", row.names = FALSE, col.names = TRUE)


