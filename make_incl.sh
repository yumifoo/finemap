trait=$1

# make .incl file
cd /fs/projects/ukbb/yu/pheno/BasicQT
/apps/statistics2/R-3.6.0/bin/Rscript --vanilla make_incl.R $trait # check directory

