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
#include "Noisy_Dynamics_Methods.hpp"

void Sample_States(double T, int Tid, std::string filename, int network_id, int world_rank, int numsteps)
{
	RegNetwork R = 	ReadTopologyFile(filename + "_" + std::to_string(network_id) + ".ids", filename + "_" + std::to_string(network_id) + ".topo", filename + "_" + std::to_string(network_id) + ".phs");
	if(!R.isvalid)
	{
		std::cout << "Network input failed. Your input files disrespected Tuco Salamanca." << "\n";
		return;
	}

	std::vector<int> state;

	std::ofstream f("outputfiles/RUN" + std::to_string(Tid) + "/frustration_output_" + std::to_string(world_rank) + ".log");
	if(!static_cast<bool>(f))
	{
		std::cout << "The output file specified is AWOL. Go figure." << "\n";
		return;
	}
	std::ofstream g("outputfiles/RUN" + std::to_string(Tid) + "/phenotype_output_" + std::to_string(world_rank) + ".log");
	if(!static_cast<bool>(g))
	{
		std::cout << "The output file specified is AWOL. Go figure." << "\n";
		return;
	}

	for(int i = 0; i < numsteps; i++)
	{
		state = GetRandomState(R.numnodes);
		Metropolis(T, 5000, state, R, &f, &g);
	}

	f.close();
	g.close();

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
	double T = std::stod(argv[3]);
	int Tid = std::stoi(argv[4]);

	Sample_States(T, Tid, "inputfiles/" + filename, network_id, world_rank, 500);

	MPI_Finalize();

	return(0);
}
