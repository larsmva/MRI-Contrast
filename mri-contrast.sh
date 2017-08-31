#!/bin/bash
# Job name:
#SBATCH --job-name=kent
#
# Project:
#SBATCH --account=nn9279k
# Wall clock limit:
#SBATCH --time='24:00:00'
#
# Max memory usage per task:
#SBATCH --mem-per-cpu=3800M
#
# Number of tasks (cores):
##SBATCH --nodes=1 --ntasks=15
#SBATCH --ntasks=2
##SBATCH --hint=compute_bound
#SBATCH --cpus-per-task=1

#SBATCH --partition=long
##SBATCH --output=output.$SCRATCH 

## Set up job environment
source /cluster/bin/jobsetup

echo $1 $2 $3

#module load gcc/4.9.2
#module load openmpi.gnu/1.8.4
#source ~oyvinev/intro/hashstack/fenics-1.5.0.abel.gnu.conf
#source ~oyvinev/fenics1.6/fenics1.6
source ~johannr/fenics-dev-2016.04.06.abel.gnu.conf

# Expand pythonpath with locally installed packages
export PYTHONPATH=$PYTHONPATH:$HOME/.local/lib/python2.7/site-packages/

# Define what to do when job is finished (or crashes)
cleanup "mkdir -p $HOME/results"
cleanup "cp -r $SCRATCH/nscm* $HOME/results" 

echo "SCRATCH is $SCRATCH"

# Copy necessary files to $SCRATCH
cp mri-contrast.py parenchyma.h5 parenchyma.xdmf "MRI0.h5" "MRI1.h5" $SCRATCH

cd $SCRATCH
ls
echo $SCRATCH
mpirun --bind-to none python nscm.py No_refinements=$1 dt_val=$2 D_val=$3

