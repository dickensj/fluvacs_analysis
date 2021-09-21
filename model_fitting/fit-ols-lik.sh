#!/bin/bash
# The interpreter used to execute the script

#“#SBATCH” directives that convey submission options:

#SBATCH --job-name=run_sims
#SBATCH --mail-user=josephdi@umich.edu
#SBATCH --mail-type=END
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=1000m 
#SBATCH --time=4:00:00
#SBATCH --account=jzelner1
#SBATCH --partition=standard
#SBATCH --output=logs/%x-%j.log


seed=124
for j in {1..12}
do
models/random_waning_ols method=sample algorithm=hmc metric=diag_e engine=nuts max_depth=10 num_warmup=1000 num_samples=1000 adapt engaged=1 delta=0.8 data file="fluvacs.data.R" random seed=$seed id=$j output refresh=25 file="results/ols_age_wane_ek_chain_${j}.csv" &

done
wait