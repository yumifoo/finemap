#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(stringr)

trait <- as.character(args[1])
chr <- 23
maf <- 0.01

source("/fs/projects/ukbb/yu/src/peak_picker_function_v2.R")
source("/fs/projects/ukbb/yu/src/merge_region_v2.R")
source("/fs/projects/ukbb/yu/src/BOLT2z.R")

##################################################
# Read in results                               #
# Adding MAF column and filtering out MAF<0.01  #
##################################################
res <- fread(sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/BOLT_sex_combined_s340951_WhiteBritish.%s.bgen.stats.gz", trait))
res <- subset(res, CHR == chr)
res$MAF <- ifelse(res$A1FREQ < 0.5, res$A1FREQ, 1-res$A1FREQ)
res <- res[res$MAF > maf, ]
res$P_BOLT_LMM_INF <- as.numeric(res$P_BOLT_LMM_INF)
dat <- res[,.(SNP, CHR, BP, CHISQ_BOLT_LMM_INF, P_BOLT_LMM_INF)]
setnames(dat, old = "CHISQ_BOLT_LMM_INF", new = "CHISQ")
setnames(dat, old = "P_BOLT_LMM_INF", new = "PVAL")

########### non-PAR #############
non.PAR.list <- fread("/fs/projects/ukbb/yu/ukb_imp_chrX_snplist.txt")
res.chr <- subset(res, BP %in% non.PAR.list$position)
dat.chr <- subset(dat, BP %in% non.PAR.list$position)

# step 1 pick significant regions
region <- peak_picker(dat.chr, window_mb = 1, plot_step = FALSE)
if(!is.null(nrow(region))){
  region <- region[,c(2,4),drop = FALSE]
}
print(paste("regions defined for", trait, "chr 23, non-PAR"))
print(region)
# step 2 merge overlapping regions
merged <- merge_region(region, window_mb = 1, data = dat.chr)
# plot the results
png(sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/plots/merge_region_sex_combined_%s_chr%d.png",
            trait, chr, trait, chr))
plot(dat.chr$BP, -log10(dat.chr$PVAL), pch = 19, col ="grey",cex=0.5, 
     xlab = "BP", ylab = "-log10(p-value)",
     main = paste("window size: 1 Mb"))
abline(h = -log10(5e-8), col = "red")
if(!is.null(nrow(merged))){
       rect(merged[,1],0,merged[,2],350,density = 30)
   rect(region[,1], -6, region[,2], -1, density = 30)
}
if(is.null(nrow(merged)) & length(merged)!= 0){
        rect(merged[1],0,merged[2],350,density = 30)
    rect(region[,1], -6, region[,2], -1, density = 30)
}
dev.off()

print("Regions merged")
print(merged)

# save the region positions 
write.table(merged, file = sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/regions/merged_region_1mb.txt",
                                   trait, chr),row.names = FALSE, col.names = FALSE, quote = FALSE)
# step 3 take the SNP with the smallest P value in each region
top.snp <- c()
if(!is.null(nrow(merged))){
  for (zz in 1:nrow(merged)) {
    stats <- subset(res.chr, BP > merged[zz,1] & BP < merged[zz,2])
    smallP <- stats[which.min(stats$P_BOLT_LMM_INF),]
    top.snp <- rbind(top.snp, smallP)
  }
  print(paste("top SNP for chr 23, non-PAR", trait))
  print(top.snp)
}

if(is.null(nrow(merged)) & length(merged)!= 0){
    stats <- subset(res.chr, BP > merged[1] & BP < merged[2])
    smallP <- stats[which.min(stats$P_BOLT_LMM_INF),]
    top.snp <- rbind(top.snp,smallP)
}
print(paste("top SNP for chr", chr, trait))
print(top.snp)


BOLT2z(top.snp, sprintf("/fs/projects/ukbb/yu/BOLT_biomarkers_agesq/%s/sex_combined/chr%d/z_files/top_snp_sex_combined_%s_chr%d_wim_1mb",
                        trait, chr, trait, chr))

########### PAR #############
#PAR.list <- fread("/fs/projects/ukbb/yu/ukb_imp_chrXY_snplist.txt")
#res.chr <- subset(res, BP %in% PAR.list$position)
#dat.chr <- subset(dat, BP %in% PAR.list$position)

# step 1 pick significant regions
#region <- peak_picker(dat.chr, window_mb = 1, plot_step =FALSE)
#if(!is.null(nrow(region))){
#  region <- region[,c(2,4),drop = FALSE]
#}
#print(paste("regions defined for", trait, "chr 23, non-PAR"))
#print(region)

# step 2 merge overlapping regions
#merged <- merge_region(region, window_mb = 1, data=dat.chr)
# plot the results
#png(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/plots/merge_region_sex_combined_%s_PAR.png",
#             trait, trait))
#plot(dat.chr$BP, -log10(dat.chr$PVAL), pch = 19, col ="grey",cex=0.5, 
#     xlab = "BP", ylab = "-log10(p-value)",
#     main = paste("window size: 1 Mb"))
#abline(h = -log10(5e-8), col = "red")
#rect(merged[,1],0,merged[,2],350,density = 30)
#rect(region[,1], -6, region[,2], -1, density = 30)
#dev.off()
#
#print("Regions merged")
#print(merged)

# save the region positions 
#write.table(merged, file = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/regions/merged_region_1mb.txt",
#                                   trait),row.names = FALSE, col.names = FALSE, quote = FALSE)

# step 3 take the SNP with the smallest P value in each region
#top.snp <- c()
#if(!is.null(nrow(merged))){
#  for (zz in 1:nrow(merged)) {
#    stats <- subset(res.chr, BP > merged[zz,1] & BP < merged[zz,2])
#    smallP <- stats[which.min(stats$P_BOLT_LMM_INF),]
#    top.snp <- rbind(top.snp, smallP)
#  }
#  print(paste("top SNP for chr 23, PAR", trait))
#  print(top.snp)
#}
#BOLT2z(top.snp, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/PAR/z_files/top_snp_sex_combined_%s_PAR_wim_1mb",
#                        trait, trait))


