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

void Sample_Mutated_Trajectories(std::string filename, int network_id, int world_rank, int numsteps)
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

	std::ofstream networkfile;

	std::vector<std::vector<int>> SS_Dec;
	std::vector<int> state, SS_F, SS_score;

	for(int mutcount = 0; mutcount < 500; mutcount++)
	{
		if(mutcount % 5 == 0)
		{
			networkfile.open("outputfiles/network_" + std::to_string(world_rank) + "_" + std::to_string(mutcount) + ".topo");
			Print_Network(R, &networkfile);
			networkfile.close();
		}

		for(int i = 0; i < numsteps; i++)
		{
			state = GetRandomState(R.numnodes);
			state = TrajectoryAsynUpdate(5000, state, R);
			SS_Dec.push_back(state);
			SS_F.push_back(CalculateFrustration(state, R));
			SS_score.push_back(Get_Phenotypic_Score(state, R));
		}
		for(int i = 0; i < SS_F.size(); i++)
		{
			f << SS_F[i] << " ";
		}
		f << "\n";
		for(int i = 0; i < SS_score.size(); i++)
		{
			g << SS_score[i] << " ";
		}
		g << "\n";
		for(int i = 0; i < SS_F.size(); i++)
		{
			for(int j = i + 1; j < SS_F.size(); j++)
			{
				h << StateOverlap(SS_Dec[i], SS_Dec[j], R) << " ";
			}
		}
		h << "\n";

		Mutate_Network(&R);

		SS_Dec.clear();
		SS_F.clear();
		SS_score.clear();
	}

	f.close();
	g.close();
	h.close();

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

	Sample_Mutated_Trajectories("inputfiles/" + filename, network_id, world_rank, 500);

	MPI_Finalize();

	return(0);
}
