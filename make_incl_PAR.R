#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

trait <- as.character(args[1])

library(data.table)

sample.XY <- fread("/fs/projects/ukbb/geno_v3/ukb22627_imp_chrXY_v3_s486429.sample")[-1,1]

x <- fread("sex_combined_basic_quantitative_traits_BOLT_s343695_agesq.txt")

trait.col <- paste0(trait,".0.0")

incl <- x[,c("IID","FID",trait.col), with = FALSE]
incl <- incl[complete.cases(incl),]

incl <- subset(incl, IID %in% sample.XY$ID_1)

write.table(incl$IID,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/sex_combined_BOLT_PAR_%s.incl",trait),
            row.names = FALSE, col.names = FALSE, quote = FALSE)