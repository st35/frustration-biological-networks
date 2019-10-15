#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include "Network_Type.hpp"
#include "Process_Topology_Methods.hpp"
#include "BFS.hpp"

void Calculate_GRC(std::string filename, int network_id, int world_rank, std::ofstream *f)
{
	RegNetwork R = 	ReadTopologyFile(filename + "_" + std::to_string(network_id) + ".ids", filename + "_" + std::to_string(network_id) + ".topo", filename + "_" + std::to_string(network_id) + ".phs");
	if(!R.isvalid)
	{
		std::cout << "Network input failed. Your input files disrespected Tuco Salamanca." << "\n";
		return;
	}

	(*f) << GRC(R) << "\n";

	return;
}

int main(int argc, char *argv[])
{
	int world_rank = 0;

	std::string filename = argv[1];

	std::ofstream f("outputfiles/GRC.log");
	for(int i = 0; i < 750; i++)
	{
		Calculate_GRC("inputfiles/" + filename, i, world_rank, &f);
	}
	f.close();

	return(0);
}
