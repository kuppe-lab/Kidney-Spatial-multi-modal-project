#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=400G
#SBATCH --array=0-14
#SBATCH --output=/beegfs/home/mgrasshoff/logs/Spatial_Clones_Kuppe/Kidney-Spatial-multi-modal-project/output.%J.%x.txt
#SBATCH --job-name=00600_MSS_Variant_Burden_Per_Spot
#SBATCH --mail-type=END
#SBATCH --mail-user=mg.slurm@proton.me
#SBATCH --time=10-00:00:00

echo "Running script:"
echo "/beegfs/data/CostaLab/Spatial_ATAC_Analysis_KuppeLab/results/SpatialClones/Kidney-Spatial-multi-modal-project/04_Mitochondrial_Clone_Identification/00600_MSS_Variant_Burden_Per_Spot.sh"
echo "Starting time:"
date

echo "Loading conda environment."
module unload R
module unload scRNA
. /beegfs/home/mgrasshoff/miniconda3/etc/profile.d/conda.sh
conda activate /beegfs/data/users/mgrasshoff/spatial_sigurd_py_env

echo "Doing array number ${SLURM_ARRAY_TASK_ID}."

echo "Loading the data."
python "/beegfs/data/CostaLab/Spatial_ATAC_Analysis_KuppeLab/results/SpatialClones/Kidney-Spatial-multi-modal-project/04_Mitochondrial_Clone_Identification/00600_MSS_Variant_Burden_Per_Spot.py" \
       ${SLURM_ARRAY_TASK_ID}
