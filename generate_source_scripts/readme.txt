# -----------------------------------------
/2020_08_03_EPI_generate_sim_files_YM_AM.sh
# -----------------------------------------

  Used to generate source code for the simulations involving Y-A (YM-AM) transitions. In effect, it 
  generates a folder structure where 1000 simulations per combination of a given epistasis scenario 
  (dominance, dominance × dominance, and coadaptation) and for both SA^Y and SA^A being involved in
  epistasis. This results in 1000 * 3 * 2 = 6000 independent simulations.

  The file calls on an executable 'seed_gen' which must have been created by compiling the files in
  ../seed_generator/ (e.g. using GCC, see compilation script). The 'seed_gen' executable must be 
  located in:
 
  $HOME/2019_04_SD_epistasis/seed_generator/seed_gen

  This can be modified in lines 39-42 of the script.

  Additionally, the script copies a file '2020_08_03_EPI_PARS_YM_AM.txt' to each simulation folder.
  This file contains standardized parameter values that are not varied in the simulations; it must
  be located in the same folder as from where this script is ran.

# -----------------------------------------
/2020_08_03_EPI_generate_sim_files_YM_AM2.sh
# -----------------------------------------

  Same as above, but with epsilon = 0 so that no epistasis occurs. This set is used to confirm that
  in absence of epistasis, turnovers occur in identical fashion for all Y-A scenarios.

# -----------------------------------------
/2020_08_03_EPI_generate_sim_files_YM_FD.sh
# -----------------------------------------

  Same, but for Y-W (YM-FD) transitions, and with either SA^Y or SA^W (not SA^A) being involved in
  epistasis. Requires an executable 'seed_gen' which can be created by compiling the files in 
  ../seed_generator_FD/, and must be located in:

  $HOME/2019_04_SD_epistasis/seed_generator_FD/seed_gen

  This can again be modified in lines 39-42 of the script.

  Additionally, the script copies a file '2020_08_03_EPI_PARS_YM_FD.txt' to each simulation folder.
  This file contains standardized parameter values that are not varied in the simulations; it must
  be located in the same folder as from where this script is ran.

# -----------------------------------------
/2020_08_03_EPI_generate_sim_files_YM_FD2.sh
# -----------------------------------------

  Same as above, but with epsilon = 0 so that no epistasis occurs. This set is used to confirm that
  in absence of epistasis, turnovers occur in identical fashion for all Y-W scenarios.

# ----------------------------
/2020_08_03_EPI_PARS_YM_AM.txt
# ----------------------------

  Contains standardized parameter values for Y-A (YM-AM) simulations. Is used in 
  /2020_08_03_EPI_generate_sim_files_YM_AM.sh

# ----------------------------
/2020_08_03_EPI_PARS_YM_FD.txt
# ----------------------------

  Contains standardized parameter values for Y-W (YM-FD) simulations. Is used in 
  /2020_08_03_EPI_generate_sim_files_YM_FD.sh