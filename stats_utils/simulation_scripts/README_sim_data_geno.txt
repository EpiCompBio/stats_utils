Regenerate molecular pheno data with e.g.:
python simulate_cont_var.py --createDF --sample-size=10000 --var-size=2000 -O cont_var_sim_data

For cov sim data do e.g.:
python simulate_cont_var.py --createDF --sample-size=10000 --var-size=10 -O cov_sim_data

For genotype data, run plink then bash conversion script:
plink --dummy 10000 1000 scalar-pheno
bash plink_to_geno.sh plink.matrixQTL plink.A-transpose plink.A-transpose.matrixQTL.geno

plink needs to be installed separately.

