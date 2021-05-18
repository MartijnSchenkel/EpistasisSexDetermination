#!/bin/bash

for n in {0..999} 
do

cat << EOT >> ${n}_FD.sh
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --partition=regular
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2000MB
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=m.a.schenkel@rug.nl

DATE=2020_08_03

NEWSD=(FD)
PREFIX=(COA DOM DXD)
INTERACT=(2 6)
LOCUS=(SAY SAII)

module load R/3.6.1-foss-2018a

for S in {0..0}
do 

SUFFIX=YM_\${NEWSD[\$S]}

for P in {0..2}
do

for L in {0..1}
do

SIM=\${PREFIX[\$P]}_\${SUFFIX}_\${LOCUS[\$L]}_${n}

cd $HOME/2019_04_SD_epistasis/\${DATE}/\${SIM}/
Rscript \${DATE}_EPI_\${SIM}.R
done
done
done
EOT
sbatch ${n}_FD.sh
rm ${n}_FD.sh
done

