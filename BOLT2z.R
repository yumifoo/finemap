BOLT2z <- function(BOLT_res, output_file_name){
  # BOLT_res is the GWAS summary statistics from BOLT-LMM
  # output_file_name is the filename of the output, should be a string of character
  
  # if using bgen support of ldstore/finemap
  #the chromosome name should be the same as in bgen file, i.e. chromosome 1 should be 01
  BOLT_res$CHR <- sprintf("%02d",BOLT_res$CHR)
  BOLT_res[BOLT_res$CHR == 23,]$ CHR <- rep("X", nrow(subset(BOLT_res, CHR == 23)))
  
  z <- data.frame(rsid = BOLT_res$SNP,
                  chromosome = BOLT_res$CHR,
                  position = BOLT_res$BP,
                  allele1 = BOLT_res$ALLELE1,
                  allele2 = BOLT_res$ALLELE0,
                  maf = BOLT_res$MAF,
                  beta = BOLT_res$BETA,
                  se = BOLT_res$SE)
  write.table(z, paste0(output_file_name,".z"), col.names = TRUE, row.names = FALSE, quote = FALSE)
}