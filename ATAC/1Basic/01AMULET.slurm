#! /bin/bash
#SBATCH --job-name=amulet
#SBATCH --nodes=1
#SBATCH --cpus-per-task 3
#SBATCH --time=24:00:00
#SBATCH -e ../tmp/amulet.err
#SBATCH -o ../tmp/amulet.out
#SBATCH --partition=basic

source /home/gaojie/miniforge3/etc/profile.d/conda.sh
conda activate scanpy
for sample in Adjacent Distant Sulcus
do
/home/gaojie/softwares/AMULET/AMULET.sh \
                                        sample_data/${sample}/fragments.tsv.gz \
                                        sample_data/${sample}/singlecell.csv \
                                        softwares/AMULET/human_autosomes.txt \
                                        softwares/AMULET/RepeatMasker_hg38_UCSC.bed \
                                        workspace/Mida_collab/AMULET/${sample} \
                                        /home/gaojie/softwares/AMULET/
done