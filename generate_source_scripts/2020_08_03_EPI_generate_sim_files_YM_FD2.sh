#!/bin/bash
# #SBATCH --nodes=1
# #SBATCH --time=1-00:00:00
# #SBATCH --partition=regular
# #SBATCH --cpus-per-task=1
# #SBATCH --mem-per-cpu=2000MB
# #SBATCH --mail-type=FAIL
# #SBATCH --mail-user=m.a.schenkel@rug.nl

module load GCC/8.2.0-2.31.1
DATE=2020_08_03
mkdir $HOME/2019_04_SD_epistasis/${DATE}

NEWSD=(FD)
PREFIX=(COA DOM DXD)
INTERACT=(2 6)
LOCUS=(SAY SAII)

for S in {0..0}
do 

SUFFIX=YM_${NEWSD[$S]}

for P in {0..2}
do

for L in {0..1}
do

for n in {1000..1999} 
do

SIM=${PREFIX[$P]}_${SUFFIX}_${LOCUS[$L]}_${n}

mkdir $HOME/2019_04_SD_epistasis/${DATE}/${SIM}

# Generate seed, what a convoluted way. If only R had a decent integer generator...
cd /home/p275703/2019_04_SD_epistasis/seed_generator_FD/
/home/p275703/2019_04_SD_epistasis/seed_generator_FD/seed_gen
cd $HOME
cp /home/p275703/2019_04_SD_epistasis/seed_generator_FD/seed.txt $HOME/2019_04_SD_epistasis/${DATE}/${SIM}/
  
FILE=${DATE}_EPI_${SIM}.R

cat << EOT >> ${FILE}
# ============#
# Information #
# ============#

# Source code accompanying the following manuscript:
# Manuscript: 	Epistatic interactions between the  sex chromosomes and autosomes can affect the 
#               stability of sex determination systems.
# Authors:    	Martijn A. Schenkel, Leo W. Beukeboom & Ido Pen.
# Affiliation:	Groningen Institute for Evolutionary Life Sciences
#               University of Groningen
#               PO Box 11103
#               9700CC Groningen, The Netherlands
# Contact:	    m.a.schenkel@rug.nl
#               i.r.pen@rug.nl

# NOTE: Must be combined with a parameter file "parameters.txt"

# Clear workspace before simulations.
rm(list=ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ================= #
# Simulation set up #
# ================= #
library(tidyverse)

# Values to be randomized. Clock used as seed generator.
seed <- as.integer(read.table("seed.txt", T))
set.seed(seed)
SAY    <- runif(1, min = 0, max = 0.05)
SAI    <- 0
SAII   <- runif(1, min = 0, max = 0.05)
EpsVal <- 0

# Define which locus interacts with EPI
# 2 = SAY
# 4 = SAI
# 6 = SAII
InteractingLocus <- ${INTERACT[$L]}

# Define file name prefix and suffix, used to separate e.g. result files from different SD transitions
Prefix <- "${PREFIX[$P]}" # generally used to denote epistasis scenario
Suffix <- "${SUFFIX}_${LOCUS[$L]}" # generally used to denote SD transition scenario


# Read parameter file and store copy in results folder
Par <- as.matrix(read.table("${DATE}_EPI_PARS_${SUFFIX}.txt", header = T))
write.table(file = "parameters.txt", x = Par)

# FileResults <- paste(${n}, ".txt", sep = "")

# ================ #
# Model parameters #
# ================ #
runtime <- Par[1] # number of generations to run the model
genskip <- Par[2] # write output every x generations
EpiEffect <- 1 + EpsVal

# recombination rates
rec <- Par[3:5]

# Frequency of focal pOne at all loci for [1] maternally and [2] paternally inherited copies
pOne <- array(data = Par[6:17], dim = c(3,2,2))
pEpi <- Par[18:19]

# Fitness and dominance parameters for all loci
sSA <- matrix(c(SAY, SAI, -SAII, -SAY, -SAI, SAII), ncol = 2, nrow = 3)
hSA <- matrix(Par[26:31], ncol = 2)

# Mutation parameters
muFD <- Par[32]
muAM <- Par[33]
muGeneration <- Par[34]

# Fitness scores per locus per sex
wXYM <- c(1, 1 + hSA[1,1] * sSA[1,1], 1 + sSA[1,1])
wIM  <- c(1, 1 + hSA[2,1] * sSA[2,1], 1 + sSA[2,1])
wIIM <- c(1, 1 + hSA[3,1] * sSA[3,1], 1 + sSA[3,1])

wXYF <- c(1, 1 + hSA[1,2] * sSA[1,2], 1 + sSA[1,2])
wIF  <- c(1, 1 + hSA[2,2] * sSA[2,2], 1 + sSA[2,2])
wIIF <- c(1, 1 + hSA[3,2] * sSA[3,2], 1 + sSA[3,2])


# ================================================ #
# Create genotype 2 gamete  transition (G2G) array  #
# ================================================ #

# Calculate proportion of all 4 haplotype produced by all 16 possible haplotype combinations.
# 1 = 00, 2 = 01, 3 = 10, 4 = 11.
# first digit is SD allele, second digit is SA allele of haplotype. 
pHaplo <- array(dim=c(3,4,4,4))
for(i in 1:length(rec))
{
  pHaplo[i,1,1,] <- c(1, 0, 0, 0)
  pHaplo[i,1,2,] <- c(1/2, 1/2, 0, 0) 
  pHaplo[i,1,3,] <- c(1/2, 0, 1/2, 0)
  pHaplo[i,1,4,] <- c((1-rec[i])/2, rec[i]/2, rec[i]/2, (1-rec[i])/2)
  pHaplo[i,2,1,] <- c(1/2, 1/2, 0, 0)
  pHaplo[i,2,2,] <- c(0, 1, 0, 0)
  pHaplo[i,2,3,] <- c(rec[i]/2, (1-rec[i])/2, (1-rec[i])/2, rec[i]/2)
  pHaplo[i,2,4,] <- c(0, 1/2, 0, 1/2)
  pHaplo[i,3,1,] <- c(1/2, 0, 1/2, 0)
  pHaplo[i,3,2,] <- c(rec[i]/2, (1-rec[i])/2, (1-rec[i])/2, rec[i]/2)
  pHaplo[i,3,3,] <- c(0, 0, 1, 0)
  pHaplo[i,3,4,] <- c(0, 0, 1/2, 1/2)
  pHaplo[i,4,1,] <- c((1-rec[i])/2, rec[i]/2, rec[i]/2, (1-rec[i])/2)
  pHaplo[i,4,2,] <- c(0, 1/2, 0, 1/2)
  pHaplo[i,4,3,] <- c(0, 0, 1/2, 1/2)
  pHaplo[i,4,4,] <- c(0, 0, 0, 1)
}

# Set frequency of EPI allele in all genotypes.
pHapEpi <- array(dim = c(2,2,2))
pHapEpi[1,1,] <- c(1, 0)
pHapEpi[1,2,] <- c(1/2, 1/2)
pHapEpi[2,1,] <- c(1/2, 1/2)
pHapEpi[2,2,] <- c(0, 1)

# Create G2G array
G2G <- array(dim=c(4,4,4,2,4,4,4,2,4,4,4,2))
for(i in 1:4)
{
  for(j in 1:4)
  {
    for(k in 1:4)
    {
      for(l in 1:2)
      {
        for(m in 1:4)
        {
          for(n in 1:4)
          {
            for(o in 1:4)
            {
              for(p in 1:2)
              {
                for(q in 1:4)
                {
                  for(r in 1:4)
                  {
                    for(s in 1:4)
                    {
                      for(t in 1:2)
                        G2G[i,j,k,l,m,n,o,p,q,r,s,t] <- 
                          pHaplo[1,i,m,q] * 
                          pHaplo[2,j,n,r] * 
                          pHaplo[3,k,o,s] * 
                          pHapEpi[l,p,t]
                    }
                  }
                }
              }   
            }
          }
        }
      }
    }
  }
}

# Function to extract gamete frequencies from G2G array.
MakeGametes <- function(G){
  F0[G[1], G[2], G[3], G[4], G[5], G[6], G[7], G[8]] * G2G[G[1], G[2], G[3], G[4], G[5], G[6], G[7], G[8],,,,]
}

# Functions to mutate gametes
MutateAM <- function(G) {
  G[,3,,] <- G[,1,,] * muAM
  G[,,1,] <- G[,,1,] * (1 - muAM)

  G[,4,,] <- G[,2,,] * muAM
  G[,,1,] <- G[,,1,] * (1 - muAM)
  G
}
MutateFD <- function(G) {
  G[,,3,] <- G[,,1,] * muFD
  G[,,1,] <- G[,,1,] * (1 - muFD)

  G[,,4,] <- G[,,2,] * muFD
  G[,,2,] <- G[,,2,] * (1 - muFD)
  G
}

# ============================================= #
# Calculate genotype frequencies for chromosome #
#     pair XY, I, and II (GS) and III (GSE)     #
# ============================================= #
GS <- array(dim = c(3,2,4))
for(i in 1:3)
{
  GS[i,1,] <- c((1 - pOne[i,1,1]) * (1 - pOne[i,2,1]), 
                (1 - pOne[i,1,1]) * (pOne[i,2,1]), 
                (pOne[i,1,1]) * (1 - pOne[i,2,1]), 
                (pOne[i,1,1]) * (pOne[i,2,1]))
  GS[i,2,] <- c((1 - pOne[i,1,2]) * (1 - pOne[i,2,2]),
                (1 - pOne[i,1,2]) * (pOne[i,2,2]), 
                (pOne[i,1,2]) * (1 - pOne[i,2,2]), 
                (pOne[i,1,2]) * (pOne[i,2,2]))
}
GSE <- array(dim = c(2,2))
GSE[1,] <- c((1 - pEpi[1]) * (1 - pEpi[2]),
             (1 - pEpi[1]) * (pEpi[2]))
GSE[2,] <- c((pEpi[1]) * (1 - pEpi[2]),
             (pEpi[1]) * (pEpi[2]))

# ======================================== #
# Initialize population and other relevant #
# arrays (fitness scores, data extraction) #
# ======================================== #

F0 <- array(dim = c(4, 4, 4, 2, 4, 4, 4, 2))
nOneM <- array(dim=c(7, 4, 4, 4, 2, 4, 4, 4, 2))
nOneP <- nOneM

nSD <- c(0,0,1,1)
nSA <- c(0,1,0,1)
nEpi <- c(0,1)

for(i in 1:4)
{
  for(j in 1:4) 
  {
    for(k in 1:4)
    {
      for(l in 1:2)
      {
        for(m in 1:4)
        {
          for(n in 1:4)
          {
            for(o in 1:4)
            {
              for(p in 1:2)
              {
                # Initialize population.
                F0[i,j,k,l,m,n,o,p] <- GS[1,1,i] * GS[2,1,j] * GS[3,1,k] * 
                  GS[1,2,m] * GS[2,2,n] * GS[3,2,o] * GSE[l,p]
                
                # Count number of focal alleles per genotype.
                nOneM[1,i,j,k,l,m,n,o,p] <- nSD[i]
                nOneM[2,i,j,k,l,m,n,o,p] <- nSA[i]
                nOneM[3,i,j,k,l,m,n,o,p] <- nSD[j]
                nOneM[4,i,j,k,l,m,n,o,p] <- nSA[j]
                nOneM[5,i,j,k,l,m,n,o,p] <- nSD[k]
                nOneM[6,i,j,k,l,m,n,o,p] <- nSA[k]
                nOneM[7,i,j,k,l,m,n,o,p] <- nEpi[l]
                
                nOneP[1,i,j,k,l,m,n,o,p] <- nSD[m]
                nOneP[2,i,j,k,l,m,n,o,p] <- nSA[m]
                nOneP[3,i,j,k,l,m,n,o,p] <- nSD[n]
                nOneP[4,i,j,k,l,m,n,o,p] <- nSA[n]
                nOneP[5,i,j,k,l,m,n,o,p] <- nSD[o]
                nOneP[6,i,j,k,l,m,n,o,p] <- nSA[o]
                nOneP[7,i,j,k,l,m,n,o,p] <- nEpi[p]
              }
            }
          }
        }
      }
    }
  }
}
nOne <- nOneM + nOneP  

# Create arrays with sex determination info
IsMale <- ((nOne[1,,,,,,,,] > 0 | nOne[3,,,,,,,,] > 0) & nOne[5,,,,,,,,] == 0)
IsFemale <- 1 - IsMale

# Definitions for different epistasis scenarios considered here.
# Dominance
EpistasisDOM <- (nOne[InteractingLocus,,,,,,,,] > 0 & nOne[7,,,,,,,,] > 0)

# Dominance x Dominance
EpistasisDXD <- (nOne[InteractingLocus,,,,,,,,] == 1 & nOne[7,,,,,,,,] == 1) 

# Co-adaptive
EpistasisCOA <- (nOne[InteractingLocus,,,,,,,,] == 0 & nOne[7,,,,,,,,] == 0) | (nOne[InteractingLocus,,,,,,,,] == 2 & nOne[7,,,,,,,,] == 2) 


# Define which type of epistasis is used. Can be changed to modify
Epistasis <- Epistasis${PREFIX[$P]}

# Calculate male and female fitness values of all genotypes.
WM <- wXYM[nOne[2,,,,,,,,] + 1] * wIM[nOne[4,,,,,,,,] + 1] * wIIM[nOne[6,,,,,,,,] + 1]
WF <- wXYF[nOne[2,,,,,,,,] + 1] * wIF[nOne[4,,,,,,,,] + 1] * wIIF[nOne[6,,,,,,,,] + 1]
W <- rep(0, length(WM))
for(i in 1:length(WF))
{
  W[i] <- max(WM[i] * IsMale[i] * (EpiEffect ^ Epistasis[i]), WF[i] * IsFemale[i])
}


Alleles <- c("YM", "SAY", "AM", "SAI", "FD", "SAII", "EPI")
LG <- c("XY", "XY", "I", "I", "II", "II", "III")
Type <- c("SD", "SA & EPI", "SD", "SA & EPI", "SD", "SA & EPI", "SA & EPI")


#=============================================================================#
#============================#                   #============================#
#============================#   Start of loop   #============================#
#============================#                   #============================#
#=============================================================================#
for(t in 2:runtime)
{
  # Which genotypes are present in the population and which of them are male and which are female?
  AreGuy <- which(F0>0 & IsMale == T,arr.ind=T,useNames=F)
  AreGirl <- which(F0>0 & IsFemale == T,arr.ind=T,useNames=F)
  
  Sperm <- rep(0, 4 * 4 * 4 * 2)
  Oocytes <- rep(0, 4 * 4 * 4 * 2)
  
  # Make gametes
  for(i in 1:nrow(AreGuy))
  {
    Sperm <- Sperm + MakeGametes(AreGuy[i,])
  }
  
  
  for(i in 1:nrow(AreGirl))
  {
    Oocytes <- Oocytes + MakeGametes(AreGirl[i,])
  }
  
  # Mutate gametes
  if(t==muGeneration)
  {
    Sperm <- MutateAM(Sperm)
    Sperm <- MutateFD(Sperm)

    Oocytes <- MutateAM(Oocytes)
    Oocytes <- MutateFD(Oocytes)
  }
  
  # Reproduce & normalize
  F1 <- Oocytes %o% Sperm
  F1 <- F1/sum(F1)
  
  
  F0 <- F1 * W
}
#=============================================================================#
#============================#                   #============================#
#============================#    End of loop    #============================#
#============================#                   #============================#
#=============================================================================#

results <- tibble(Seed = numeric(),
                  Generation = numeric(),
                  SexRatio = numeric(),
                  Frequency = numeric(),
                  Sex = character(),
                  LinkageGroup = character(),
                  Copy = character(),
                  Allele = character(),
                  Type = character(),
                  SAY = numeric(),
                  SAI = numeric(),
                  SAII = numeric(),
                  Epsilon = numeric(),
                  Scenario = character(),
                  Epistasis = character())
for(i in 1:7)
{
  temp <- tibble(Seed = seed,
                 Generation = 1,
                 SexRatio = F0 %*% IsMale,
                 Frequency = nOneM[i,,,,,,,,] %*% (F0 * IsMale),
                 Sex = "M",
                 LinkageGroup = LG[i],
                 Copy = "Maternal",
                 Allele = Alleles[i],
                 Type = Type[i],
		 SAY = SAY,
                 SAI = SAI,
                 SAII = SAII,
                 Epsilon = EpsVal,
                 Scenario = Suffix,
                 Epistasis = Prefix)
  
  results <- results %>% bind_rows(temp)
  
  temp <- tibble(Seed = seed,
                 Generation = 1,
                 SexRatio = F0 %*% IsMale,
                 Frequency = nOneM[i,,,,,,,,] %*% (F0 * IsFemale),
                 Sex = "F",
                 LinkageGroup = LG[i],
                 Copy = "Maternal",
                 Allele = Alleles[i],
                 Type = Type[i],
		 SAY = SAY,
                 SAI = SAI,
                 SAII = SAII,
                 Epsilon = EpsVal,
                 Scenario = Suffix,
                 Epistasis = Prefix)
  
  results <- results %>% bind_rows(temp)
  
  
  
  temp <- tibble(Seed = seed,
                 Generation = 1,
                 SexRatio = F0 %*% IsMale,
                 Frequency = nOneP[i,,,,,,,,] %*% (F0 * IsMale),
                 Sex = "M",
                 LinkageGroup = LG[i],
                 Copy = "Paternal",
                 Allele = Alleles[i],
                 Type = Type[i],
		 SAY = SAY,
                 SAI = SAI,
                 SAII = SAII,
                 Epsilon = EpsVal,
                 Scenario = Suffix,
                 Epistasis = Prefix)
  
  results <- results %>% bind_rows(temp)
  
  temp <- tibble(Seed = seed,
                 Generation = 1,
                 SexRatio = F0 %*% IsMale,
                 Frequency = nOneP[i,,,,,,,,] %*% (F0 * IsFemale),
                 Sex = "F",
                 LinkageGroup = LG[i],
                 Copy = "Paternal",
                 Allele = Alleles[i],
                 Type = Type[i],
		 SAY = SAY,
                 SAI = SAI,
                 SAII = SAII,
                 Epsilon = EpsVal,
                 Scenario = Suffix,
                 Epistasis = Prefix)
  
  results <- results %>% bind_rows(temp)  
}

# Save results for further analysis
write.table(file = paste(${n}, ".txt", sep = ""), x = results, row.names = F, sep = "\t")   
EOT

cp ${FILE} $HOME/2019_04_SD_epistasis/${DATE}/${SIM}
cp ${DATE}_EPI_PARS_${SUFFIX}.txt $HOME/2019_04_SD_epistasis/${DATE}/${SIM}

rm ${FILE}

done
done
done
done
