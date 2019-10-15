double CalculateEnergy(std::vector<int> state, RegNetwork R)
{
	double energy = 0.0;
	int J_ij, s_i, s_j;
	for(int i = 0; i < R.numedges; i++)
	{
		s_i = state[R.interactions[i][0]];
		s_j = state[R.interactions[i][1]];
		J_ij = R.interactions[i][2];

		if(s_i == 0)
		{
			s_i = -1;
		}
		if(s_j == 0)
		{
			s_j = -1;
		}
		if(J_ij == 2)
		{
			J_ij = -1;
		}

		energy += J_ij*s_i*s_j;
	}

	return(-energy);
}

void Metropolis(double T, int numsteps, std::vector<int> state, RegNetwork R, std::ofstream *f, std::ofstream *g)
{
	if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return;
	}

	std::vector<int> newstate;
	for(int i = 0; i < state.size(); i++)
	{
		newstate.push_back(state[i]);
	}

	int nodetoupdate, flag = 0;
	double DeltaE;
	for(int t = 0; t < numsteps; t++)
	{
		for(int i = 0; i < R.numnodes; i++)
		{
			nodetoupdate = (int) (distribution(generator)*R.numnodes);
			for(int j = 0; j < R.numnodes; j++)
			{
				if(j == nodetoupdate)
				{
					newstate[j] = (int) (!state[j]);
				}
				else
				{
					newstate[j] = state[j];
				}
			}
			DeltaE = CalculateEnergy(newstate, R) - CalculateEnergy(state, R);
			flag = 0;
			if(DeltaE < 0.0)
			{
				flag = 1;
			}
			else
			{
				if(distribution(generator) < std::exp(-DeltaE / T))
				{
					flag = 1;
				}
			}
			if(flag == 1)
			{
				state[nodetoupdate] = newstate[nodetoupdate];
			}
		}

		if(t > numsteps / 2 && t % 50 == 0)
		{
			(*f) << CalculateFrustration(state, R) << " ";
			(*g) << Get_Phenotypic_Score(state, R) << " ";
		}
	}
	(*f) << "\n";
	(*g) << "\n";

	return;
}

std::vector<int> Noisy_Asyn_Dynamics(double eta, int numsteps, std::vector<int> state, RegNetwork R, std::ofstream *f, std::ofstream *g)
{
	if(state.size() != R.numnodes)
	{
		std::cout << "The initial state is so wrong that you should stand there in your wrongness and be wrong." << "\n";
		return(std::vector<int>());
	}

	int errflag = 0;

	int nodetoupdate, input;

	int source, target, interactiontype;
	for(int runcount = 0; runcount < numsteps; runcount++)
	{
		for(int step = 0; step < R.numnodes; step++)
		{
			nodetoupdate = (int) (distribution(generator)*R.numnodes);
			input = 0;
			for(int i = 0; i < R.numedges; i++)
			{
				source = R.interactions[i][0];
				target = R.interactions[i][1];
				if(target != nodetoupdate)
				{
					continue;
				}
				interactiontype = R.interactions[i][2];
				if(interactiontype == 1)
				{
					if(state[source] == 0)
					{
						input += (1)*(-1);
					}
					else if(state[source] == 1)
					{
						input += (1)*(1);
					}
					else
					{
						errflag = 1;
					}
				}
				else if(interactiontype == 2)
				{
					if(state[source] == 0)
					{
						input += (-1)*(-1);
					}
					else if(state[source] == 1)
					{
						input += (-1)*(1);
					}
					else
					{
						errflag = 1;
					}
				}
				else
				{
					errflag = 1;
				}
			}
			if(errflag == 1)
			{
				std::cout << "Wicked weird stuff going on here. Take some time and think about what you have done." << "\n";
				return(std::vector<int>());
			}
			if(input > 0)
			{
				state[nodetoupdate] = 1;
			}
			else if(input < 0)
			{
				state[nodetoupdate] = 0;
			}
			else
			{
				state[nodetoupdate] = state[nodetoupdate];
			}

			if(distribution(generator) < eta)
			{
				state[nodetoupdate] = (int) (!state[nodetoupdate]);
			}
		}

		if(runcount > numsteps / 2 && runcount % 50 == 0)
		{
			(*f) << CalculateFrustration(state, R) << " ";
			(*g) << Get_Phenotypic_Score(state, R) << " ";
		}
	}
	(*f) << "\n";
	(*g) << "\n";

	return(state);
}
