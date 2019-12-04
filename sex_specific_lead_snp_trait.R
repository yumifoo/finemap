library(data.table)

args = commandArgs(trailingOnly=TRUE)
trait <- as.character(args[1])

chromosomes <- 1:23

female <- c()
male <- c()
combine <- c()

res.f <- fread(sprintf("BOLT_female_s184583_WhiteBritish.%s.bgen.stats.gz", trait))
res.m <- fread(sprintf("BOLT_male_s159112_WhiteBritish.%s.bgen.stats.gz",trait))
res.c <- fread(sprintf("BOLT_sex_combined_s343695_WhiteBritish.%s.bgen.stats.gz", trait))


for(ii in seq_along(chromosomes)){
  chr <- chromosomes[ii]
  lead.snp <- read.table(sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/lead_snp_cond/lead_snps_cond_chr%d.txt", trait, chr, chr),
                         header = FALSE, stringsAsFactors = FALSE)
  lead.snp.female <- subset(res.f, SNP %in% lead.snp[,1])
  write.table(lead.snp.female, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/lead_snp_cond/BOLT_female_res_lead_snps_cond_chr%d.txt",trait, chr, chr),
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  female <- rbind(female, lead.snp.female)
  
  lead.snp.male <- subset(res.m, SNP %in% lead.snp[,1])
  write.table(lead.snp.male, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/lead_snp_cond/BOLT_male_res_lead_snps_cond_chr%d.txt", trait,chr, chr),
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  male <- rbind(male, lead.snp.male)
  
  lead.snp.combine <- subset(res.c, SNP %in% lead.snp[,1])
  write.table(lead.snp.combine, sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/chr%d/lead_snp_cond/BOLT_combine_res_lead_snps_cond_chr%d.txt",trait,chr, chr),
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  combine <- rbind(combine, lead.snp.combine)
}
write.table(female,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/BOLT_female_res_lead_snps_cond.txt", trait),
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(male,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/BOLT_male_res_lead_snps_cond.txt", trait),
            row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(combine,sprintf("/fs/projects/ukbb/yu/BOLT_basicQT_agesq/%s/sex_combined/BOLT_combine_res_lead_snps_cond.txt", trait),
            row.names = FALSE, col.names = TRUE, quote = FALSE)
