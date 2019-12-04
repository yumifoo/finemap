trait=$1

mkdir -p $trait/{sex_combined,female,male,plots}/chr{1..23}/{regions,plots,master,z_files,ld,finemap,lead_snp_cond}
mkdir -p $trait/{sex_combined,female,male}/PAR/{regions,plots,master,z_files,ld,finemap,lead_snp_cond}

# make .incl file
cd /fs/projects/ukbb/yu/pheno/BasicQT
/apps/statistics2/R-3.6.0/bin/Rscript --vanilla make_incl.R $trait # check directory


cd /fs/projects/ukbb/yu/BOLT_basicQT_agesq
# define significant regions
# pick the top SNPs of each region to calculate ld in ldstore
for i in {1..23}
do
    /apps/statistics2/R-3.6.0/bin/Rscript --vanilla trait_region.R $trait $i
done

/apps/statistics2/R-3.6.0/bin/Rscript --vanilla trait_region_X.R $trait


/apps/statistics2/R-3.6.0/bin/Rscript --vanilla top_snp_make_master.R $trait

cd /fs/projects/ukbb/yu/BOLT_basicQT_agesq/$trait/sex_combined
# run ldstore on top SNPs
/fs/projects/ukbb/christian/binaries/150419/ldstore_v2.0b_x86_64 --in-files chr1-23_master --write-bcor --write-bdose --n-threads 8 --cpu-mem 100
/fs/projects/ukbb/christian/binaries/150419/ldstore_v2.0b_x86_64 --in-files chr1-23_master --bcor-to-text

cd /fs/projects/ukbb/yu/BOLT_basicQT_agesq
# check the results of ld store if is needed for further merge
/apps/statistics2/R-3.6.0/bin/Rscript --vanilla check_top_ld.R $trait

# make master files
/apps/statistics2/R-3.6.0/bin/Rscript --vanilla make_master.R $trait

