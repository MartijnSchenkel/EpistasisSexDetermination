#include "randomnumbers.h"
#include <fstream>
#include <iostream>
int main()
{	
	// SET UP SIMULATION //
	long seed = randomize();	
	std::ofstream sf("seed.txt");
	sf << "seed" << std::endl;
	sf << seed << std::endl;
	return(0);
}
