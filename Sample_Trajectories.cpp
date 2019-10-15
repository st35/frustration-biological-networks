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

void Sample_Trajectories(std::string filename, int network_id, int world_rank, int numsteps)
{
	RegNetwork R = 	ReadTopologyFile(filename + "_" + std::to_string(network_id) + ".ids", filename + "_" + std::to_string(network_id) + ".topo", filename + "_" + std::to_string(network_id) + ".phs");
	if(!R.isvalid)
	{
		std::cout << "Network input failed. Your input files disrespected Tuco Salamanca." << "\n";
		return;
	}

	std::vector<std::vector<int>> SS_Dec;
	std::vector<int> state, SS_F, SS_score;
	std::vector<double> SS_Q;

	for(int i = 0; i < numsteps; i++)
	{
		state = GetRandomState(R.numnodes);
		state = TrajectoryAsynUpdate(5000, state, R);
		SS_Dec.push_back(state);
		SS_F.push_back(CalculateFrustration(state, R));
		SS_score.push_back(Get_Phenotypic_Score(state, R));
	}

	std::ofstream f("outputfiles/output_" + std::to_string(world_rank) + ".log");
	if(!static_cast<bool>(f))
	{
		std::cout << "The output file specified is AWOL. Go figure." << "\n";
		return;
	}
	for(int i = 0; i < SS_F.size(); i++)
	{
		f << SS_F[i] << " " << SS_score[i] << "\n";
	}
	f.close();

	f.open("outputfiles/qab_output_" + std::to_string(world_rank) + ".log");
	if(!static_cast<bool>(f))
	{
		std::cout << "The output file specified is AWOL. Go figure." << "\n";
		return;
	}
	for(int i = 0; i < SS_F.size(); i++)
	{
		for(int j = i + 1; j < SS_F.size(); j++)
		{
			f << StateOverlap(SS_Dec[i], SS_Dec[j], R) << "\n";
		}
	}
	f.close();

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

	Sample_Trajectories("inputfiles/" + filename, network_id, world_rank, 500);

	MPI_Finalize();

	return(0);
}
