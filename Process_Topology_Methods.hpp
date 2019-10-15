RegNetwork ReadTopologyFile(std::string filename1, std::string filename2, std::string filename3)
{
	RegNetwork R;
	R.numnodes = 0;
	R.numedges = 0;
	R.nodeIDs = std::vector<std::string>();
	R.interactions = std::vector<std::vector<int>>();
	R.phenotype_pos = std::vector<int>();
	R.phenotype_neg = std::vector<int>();

	std::string nodeName;
	int nodeID, numlines;

	std::ifstream f(filename1);
	if(!static_cast<bool>(f))
	{
		std::cout << "Error reading the .ids file." << "\n";
		return(R);
	}
	while(f >> nodeName >> nodeID)
	{
		R.nodeIDs.push_back(nodeName);
		R.numnodes += 1;
	}
	f.close();

	f.open(filename1);
	numlines = 0;
	while(getline(f, nodeName))
	{
		numlines += 1;
	}
	f.close();

	if(numlines != R.nodeIDs.size())
	{
		std::cout << "The format of the .ids file is ridonc, and frankly, completely wrong." << "\n";
		return(R);
	}

	std::string source, target;
	int sourceID, targetID, interactionType;

	f.open(filename2);
	if(!static_cast<bool>(f))
	{
		std::cout << "Error reading the .topo file." << "\n";
		return(R);
	}
	getline(f, source);
	while(f >> source >> target >> interactionType)
	{
		if(std::find(R.nodeIDs.begin(), R.nodeIDs.end(), source) == R.nodeIDs.end())
		{
			std::cout << "One of the network nodes is missing from the .ids file." << "\n";
			return(R);
		}
		if(std::find(R.nodeIDs.begin(), R.nodeIDs.end(), target) == R.nodeIDs.end())
		{
			std::cout << "One of the network nodes is missing from the .ids file." << "\n";
			return(R);
		}
		sourceID = std::find(R.nodeIDs.begin(), R.nodeIDs.end(), source) - R.nodeIDs.begin();
		targetID = std::find(R.nodeIDs.begin(), R.nodeIDs.end(), target) - R.nodeIDs.begin();
		R.interactions.push_back({0, 0, 0});
		R.interactions[R.numedges][0] = sourceID;
		R.interactions[R.numedges][1] = targetID;
		R.interactions[R.numedges][2] = interactionType;
		R.numedges += 1;
	}
	f.close();

	f.open(filename2);
	numlines = 0;
	while(getline(f, source))
	{
		numlines += 1;
	}
	f.close();

	if(numlines - 1 != R.interactions.size())
	{
		std::cout << "The format of the .topo file is ridonc, and frankly, completely wrong." << "\n";
		return(R);
	}

	int nodeid_pos, nodeid_neg;
	f.open(filename3);
	if(!static_cast<bool>(f))
	{
		std::cout << "Phenotype file could not be read for this network. Moving on." << "\n";
	}
	else
	{
		while(f >> nodeid_pos >> nodeid_neg)
		{
			R.phenotype_pos.push_back(nodeid_pos);
			R.phenotype_neg.push_back(nodeid_neg);
		}
	}
	f.close();

	if(R.numnodes == R.nodeIDs.size() && R.numedges == R.interactions.size())
	{
		R.isvalid = true;
	}

	return(R);
}

void Print_Network(RegNetwork R, std::ofstream *f)
{
	(*f) << "Source\tTarget\tType\n";
	for(int i = 0; i < R.interactions.size(); i++)
	{
		(*f) << R.nodeIDs[R.interactions[i][0]] << "\t";
		(*f) << R.nodeIDs[R.interactions[i][1]] << "\t";
		(*f) << R.interactions[i][2] << "\n";
	}

	return;
}