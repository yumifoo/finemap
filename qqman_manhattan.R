#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

trait <- as.character(args[1])

library(qqman)
library(data.table)

res <- fread(sprintf("BOLT_sex_combined_s343695_WhiteBritish.%s.bgen.stats.gz", trait))

df <- data.frame(SNP = res$SNP,
                 CHR = res$CHR,
                 BP = res$BP,
                 P = as.numeric(res$P_BOLT_LMM_INF))

df[df$P<1e-323,]$P <- 1e-323

png(filename = sprintf("%s/plots/manhattan_plot_%s_sex_combined.png", trait, trait),
        width = 12, height = 7, units="in", res = 300)
manhattan(df, main = paste(trait,"sex-combined"),
                    chrlabs =c(1:22,"X"), suggestiveline = FALSE)
dev.off()
