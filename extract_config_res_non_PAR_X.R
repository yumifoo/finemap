#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
trait <- as.character(args[1])

chr <- 23

all.causal <- c()
region.file <- sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/regions/merged_region_1mb.txt",trait, chr)
if(file.exists(region.file)){
  n.region <- as.integer(system2("wc",
                                 args = c("-l",
                                          region.file,
                                          " | awk '{print $1}'"),
                                 stdout = TRUE))
  if(n.region > 0){
    file.ls <- list.files(path = sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap", trait, chr),
                          pattern = "config$")
    causal <- c()
    for(jj in seq_along(file.ls)){
      config <- read.table(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/finemap/%s",trait, chr,file.ls[jj]), header = TRUE, stringsAsFactors = FALSE)
      res <- tail(config,1)
      causal <- rbind(causal,res)
      print(paste0(jj,"/",length(file.ls)))
    }
    causal <- cbind(causal,rep(paste0("chr",chr)))
    all.causal <- rbind(all.causal, causal)
    print(paste("chr",chr))
    colnames(all.causal) <- c("rank","config","loglik","chisq","pvalue", "chromosome")
    write.table(all.causal,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/lead_snp_cond/causal_SNP_chr%d_cond_FINEMAP.txt",trait, chr, chr), row.names= FALSE, col.names = TRUE, quote = FALSE)
  }
}

