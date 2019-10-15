class RegNetwork
{
	public:
		bool isvalid = false;
		int numnodes = 0, numedges = 0;
		std::vector<std::string> nodeIDs;
		std::vector<std::vector<int>> interactions;
		std::vector<int> phenotype_pos, phenotype_neg;
};