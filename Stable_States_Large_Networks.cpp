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

std::string State_To_String(std::vector<int> state)
{
	std::string s = "";
	for(int i = 0; i < state.size(); i++)
	{
		s = s + std::to_string(state[i]);
	}

	return(s);
}

void Get_Stable_States(std::string filename, int network_id, int world_rank)
{
	RegNetwork R = 	ReadTopologyFile(filename + "_" + std::to_string(network_id) + ".ids", filename + "_" + std::to_string(network_id) + ".topo", filename + "_" + std::to_string(network_id) + ".phs");
	if(!R.isvalid)
	{
		std::cout << "Network input failed. Your input files disrespected Tuco Salamanca." << "\n";
		return;
	}

	long firststate, laststate, singlerunsize;
	int degenerate_flag = 0;
	singlerunsize = ((long) std::pow(2, 20)) / (16*16);
	if(singlerunsize == 0)
	{
		degenerate_flag = 1;
		singlerunsize = (long) std::pow(2, R.numnodes) - 1;
	}

	std::vector<int> state;
	boost::unordered_map<std::string, long> SS_Dec;
	boost::unordered_map<std::string, int> SS_F, SS_score;

	for(int i = 0; i < 16; i++)
	{	

		firststate = 16*world_rank*singlerunsize + i*singlerunsize;
		laststate = 16*world_rank*singlerunsize + i*singlerunsize + singlerunsize - 1;
		for(long j = firststate; j <= laststate; j++)
		{
			if(j > singlerunsize && degenerate_flag == 1)
			{
				break;
			}
			state = GetRandomState(R.numnodes);
			state = TrajectoryAsynUpdate(R.numnodes*R.numnodes, state, R);
			if(IsStableState(state, R))
			{
				if(SS_Dec.find(State_To_String(state)) == SS_Dec.end())
				{
					SS_Dec[State_To_String(state)] = 1;
					SS_F[State_To_String(state)] = CalculateFrustration(state, R);
					SS_score[State_To_String(state)] = Get_Phenotypic_Score(state, R);
				}
				else
				{
					SS_Dec[State_To_String(state)] += 1;
				}
			}
		}
	}

	std::ofstream f("outputfiles/RUN" + std::to_string(network_id) + "/output_" + std::to_string(world_rank) + ".log");
	if(!static_cast<bool>(f))
	{
		std::cout << "The output file specified is AWOL. Go figure." << "\n";
		return;
	}
	for(auto p: SS_Dec)
	{
		f << p.second << " " << SS_F[p.first] << " " << SS_score[p.first] << " ";
		for(int j = 0; j < p.first.length(); j++)
		{
			f << p.first[j] << " ";
		}
		f << "\n";
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

	Get_Stable_States("inputfiles/" + filename, network_id, world_rank);

	MPI_Finalize();

	return(0);
}
