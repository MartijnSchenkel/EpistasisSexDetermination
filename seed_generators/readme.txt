/seed_generator/

Contains the source code for a random number generator that is used to generate unique seeds for
each independent simulation. Source code consists of main.cpp, randomnumbers.cpp, and 
randomnumbers.h. When compiled, it generates the executable seed_gen, which can be called to 
generate a unique seed contained in the file seed.txt.

/seed_generator_FD/

Same as /seed_generator, but simply a separate instance so that seeds for Y-W (YM-FD) transition can
be generated at the same time as Y-A (YM-AM) transitions.

/compile_seed_generator.txt

Compiles the seed generator using the files in /seed_generator

/compile_seed_generator_FD.txt

Compiles the seed generator using the files in /seed_generator_FD