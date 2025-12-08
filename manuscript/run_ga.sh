#!/bin/bash
sbatch <<EOT
#!/bin/bash
#SBATCH -J ga_saem
#SBATCH -o /home/stexocae/li_lab/saem/manuscript/out/slurm/ga.o%j
#SBATCH -e /home/stexocae/li_lab/saem/manuscript/out/slurm/ga.e%j
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3

. "/home/stexocae/miniconda3/etc/profile.d/conda.sh"
conda activate saem
python3 /home/stexocae/li_lab/saem/manuscript/ga_params.py

EOT
