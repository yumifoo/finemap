trait=$1

# make .incl file
cd /fs/projects/ukbb/yu/pheno/Biomarkers
/apps/statistics2/R-3.6.0/bin/Rscript --vanilla make_incl_biomarkers.R $trait # check directory

