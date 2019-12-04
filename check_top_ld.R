#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

trait <- as.character(args[1])

chromosomes <- 1:23

need_merge <- c()

for(ii in 1:length(chromosomes)){
  chr <- chromosomes[ii]
  ld <- read.table( sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/ld/top_snp_sex_combined_win_1mb.ld",trait, chr), header = FALSE)
  ld <- ld^2
  ld <- as.matrix(ld)
  z <- read.table(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/z_files/top_snp_sex_combined_%s_chr%d_wim_1mb.z",trait, chr, trait, chr), header = TRUE, stringsAsFactors = FALSE)
  rownames(ld) <- z$rsid
  colnames(ld) <- z$rsid
  get_upper_tri<-function(cormat){
    cormat[lower.tri(cormat)] <- NA
    return(cormat)
  }
  ld_upper <- get_upper_tri(ld) 
  for (jj in 1:nrow(ld_upper)){
    ld_upper[jj,jj] <- NA
  }
  print(paste("The largest ld^2 value on chromosome", chr, "is", max(ld_upper, na.rm = TRUE)))
  print("Any LD^2 > 0.1?")
  print(any(ld_upper > 0.1, na.rm = TRUE))
  need_merge[ii] <- any(ld_upper > 0.1, na.rm = TRUE)
}

names(need_merge) <- 1:23
need_merge <- names(need_merge)[need_merge]
need_merge

for(ii in seq_along(need_merge)){
  print(paste("chr",need_merge[ii]))
  ld <- read.table(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%s/ld/top_snp_sex_combined_win_1mb.ld",trait, need_merge[ii]), header = FALSE)
  ld <- ld^2
  index <- which(ld>0.1 & ld < 1, arr.ind = TRUE)[1,]
  # THE INDEX IS THE ROW NUMBER IN THE REGION FILE

  region <- read.table(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%s/regions/merged_region_1mb.txt",trait, need_merge[ii]),header = FALSE)
  print(region)
  region <- rbind(region,c(min(region[index,]),max(region[index,])))
  print(region)
  region <- region[-index,]
  region <- region[order(region[,1]),]
  print(region)
  write.table(region, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%s/regions/merged_region_1mb.txt",trait, need_merge[ii]), quote = FALSE, row.names = FALSE, col.names = FALSE)
}
