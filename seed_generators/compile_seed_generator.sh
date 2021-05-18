#!/bin/bash
module load GCC/8.2.0-2.31.1

cd /home/p275703/2019_04_SD_epistasis/seed_generator/

g++ -x c++ /home/p275703/2019_04_SD_epistasis/seed_generator/*.cpp -x c++-header /home/p275703/2019_04_SD_epistasis/seed_generator/*.h -o /home/p275703/2019_04_SD_epistasis/seed_generator/seed_gen -std=c++17
