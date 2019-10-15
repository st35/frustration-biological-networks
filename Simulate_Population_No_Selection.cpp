#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <boost/unordered_map.hpp>
#include <mpi.h>
#include "Network_Type.hpp"
#include "Process_Topology_Methods.hpp"
#include "Boolean_Modeling_Methods.hpp"
#include "Mutate_Network.hpp"
#include "Population_Methods.hpp"

void Simulate_Population(std::string filename, int network_id, int world_rank, int numsteps)
{
	RegNetwork R = 	ReadTopologyFile(filename + "_" + std::to_string(network_id) + ".ids", filename + "_" + std::to_string(network_id) + ".topo", filename + "_" + std::to_string(network_id) + ".phs");
	if(!R.isvalid)
	{
		std::cout << "Network input failed. Your input files disrespected Tuco Salamanca." << "\n";
		return;
	}

	std::ofstream f("outputfiles/frustration_output_" + std::to_string(world_rank) + ".log");
	if(!static_cast<bool>(f))
        {
                std::cout << "The output file specified is AWOL. Go figure." << "\n";
                return;
        }
	std::ofstream g("outputfiles/phenotype_output_" + std::to_string(world_rank) + ".log");
	if(!static_cast<bool>(g))
        {
                std::cout << "The output file specified is AWOL. Go figure." << "\n";
                return;
        }
	std::ofstream h("outputfiles/qab_output_" + std::to_string(world_rank) + ".log");
	if(!static_cast<bool>(h))
        {
                std::cout << "The output file specified is AWOL. Go figure." << "\n";
                return;
        }
	std::ofstream meanfile("outputfiles/meanf_output_" + std::to_string(world_rank) + ".log");
	if(!static_cast<bool>(meanfile))
        {
                std::cout << "The output file specified is AWOL. Go figure." << "\n";
                return;
        }

	std::vector<RegNetwork> P = Generate_Same_Network_Population(R, 500);
	Simulate_No_Selection(P, 1000, 50, &f, &g, &h, &meanfile, world_rank);

	f.close();
	g.close();
	h.close();
	meanfile.close();

	return;
}

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);
    	int world_size;
    	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    	int world_rank = 0;
    	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	InitializeRandomNumGen(std::time(NULL) + world_rank);

	std::string filename = argv[1];
	int network_id = std::stoi(argv[2]);

	Simulate_Population("inputfiles/" + filename, network_id, world_rank, 5);

	MPI_Finalize();

	return(0);
}
