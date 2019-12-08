#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
trait <- as.character(args[1])

chromosomes <- 1:23

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
      cond.res <- read.table(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/lead_snp_cond/causal_SNP_chr%d_cond_FINEMAP.txt",trait, chr, chr),
                             header = TRUE, stringsAsFactors = FALSE)
      cond.res.sep <- tidyr::separate_rows(cond.res,config,sep = ",")
      lead.snp <- cond.res.sep$config
      write.table(lead.snp,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/lead_snp_cond/lead_snps_cond_chr%d.txt",trait, chr, chr), row.names= FALSE, col.names = FALSE, quote = FALSE)
    }
  }
}