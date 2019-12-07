#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

trait <- as.character(args[1])

devtools::install_github('anastasia-lucas/hudson')
library(hudson)
library(data.table)

female <- fread(sprintf("BOLT_female_s184583_WhiteBritish.%s.bgen.stats.gz",trait), colClasses = c(CHR = "character"))
male <- fread(sprintf("BOLT_male_s159112_WhiteBritish.%s.bgen.stats.gz",trait), colClasses = c(CHR = "character"))


# change 23 to X
female[female$CHR == "23",]$CHR <- rep("X", nrow(subset(female, CHR == "23")))
male[male$CHR == "23",]$CHR <- rep("X", nrow(subset(male, CHR == "23")))

gwas.f <- data.frame(SNP = female$SNP,
                     CHR = female$CHR,
                     POS = female$BP,
                     pvalue = as.numeric(female$P_BOLT_LMM_INF))
gwas.f[gwas.f$pvalue<1e-323,]$pvalue <- 1e-323

gwas.m <- data.frame(SNP = male$SNP,
                     CHR = male$CHR,
                     POS = male$BP,
                     pvalue = as.numeric(male$P_BOLT_LMM_INF))
gwas.m[gwas.m$pvalue<1e-323,]$pvalue <- 1e-323

p <- gmirror(top = gwas.f, bottom = gwas.m,
        tline = 5e-8, bline = 5e-8,
        yaxis = c(expression(paste("-log"[10], "(p-female)", sep="")), expression(paste("-log"[10], "(p-male)", sep=""))),
        opacity = 0.7,
        toptitle = paste0(trait, ", female"), bottomtitle = paste0(trait, ", male"),
        file = sprintf("%s/plots/miami_plot_%s_female_male", trait, trait)
        )
