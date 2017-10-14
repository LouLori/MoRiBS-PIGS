#!/bin/bash
#SBATCH --job-name=h2o-h2
#SBATCH --output=/work/lduranca/bosonic-system/test-run-PIMC/pimc.out
#SBATCH --time=10-00:00
#SBATCH --mem-per-cpu=1200mb
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=1
#
cp /home/lduranca/MoRiBS-PIMC/examples/bosonic-system/pes.tab /work/lduranca/bosonic-system/test-run-PIMC/
cp /home/lduranca/MoRiBS-PIMC/examples/bosonic-system/qmc.input /work/lduranca/bosonic-system/test-run-PIMC/
cp /home/lduranca/MoRiBS-PIMC/examples/bosonic-system/isoH2H208.pot /work/lduranca/bosonic-system/test-run-PIMC/
cp /home/lduranca/MoRiBS-PIMC/examples/bosonic-system/H2O* /work/lduranca/bosonic-system/test-run-PIMC/
cp /home/lduranca/MoRiBS-PIMC/pimc /work/lduranca/bosonic-system/test-run-PIMC/
#
rm -rf /scratch/lduranca/bosonic-system/test-run
mkdir -p /scratch/lduranca/bosonic-system/test-run/results
cp /work/lduranca/bosonic-system/test-run-PIMC/pes.tab /scratch/lduranca/bosonic-system/test-run
cp /work/lduranca/bosonic-system/test-run-PIMC/qmc.input /scratch/lduranca/bosonic-system/test-run
cp /work/lduranca/bosonic-system/test-run-PIMC/isoH2H208.pot /scratch/lduranca/bosonic-system/test-run
cp /work/lduranca/bosonic-system/test-run-PIMC/H2O* /scratch/lduranca/bosonic-system/test-run
cp /work/lduranca/bosonic-system/test-run-PIMC/pimc /scratch/lduranca/bosonic-system/test-run
cd /scratch/lduranca/bosonic-system/test-run
#####valgrind --leak-check=full -v --show-leak-kinds=all ./pimc
./pimc
rm -rf /work/lduranca/bosonic-system/test-run
mv /scratch/lduranca/bosonic-system/test-run /work/lduranca/bosonic-system
