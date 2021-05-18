#------------------------
/2020_08_03_run_EPI_AM.sh
#------------------------

  Generates 1000 jobscripts to run the Y-A (YM-AM) transitions. Each jobscript loops over the 
  combinations of epistasis scenario (dominance, dominance × dominance, coadaptation) and which SA
  locus (SA^Y or SA^A) is involved. To execute this script, the simulation source codes must be
  generated using the files in ../generate_source_scripts/ (for details see 'readme.txt' in this 
  folder).

#------------------------
/2020_08_03_run_EPI_AM2.sh
#------------------------

  Same as above, but for the set of Y-A simulations where epsilon = 0.

#------------------------
/2020_08_03_run_EPI_FD.sh
#------------------------

  Generates 1000 jobscripts to run the Y-W (YM-FD) transitions. Each jobscript loops over the 
  combinations of epistasis scenario (dominance, dominance × dominance, coadaptation) and which SA
  locus (SA^Y or SA^A) is involved. To execute this script, the simulation source codes must be
  generated using the files in ../generate_source_scripts/ (for details see 'readme.txt' in this 
  folder).

#------------------------
/2020_08_03_run_EPI_FD2.sh
#------------------------

  Same as above, but for the set of Y-A simulations where epsilon = 0.
