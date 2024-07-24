#include "particle.hpp"
#include "print.hpp"
#include "reaction_info.hpp"
#include "reaction_network.hpp"
#include "string_utility.hpp"

int
main()
{
	std::vector<int>             a{ 1, 2, 3, 4, 5 };
	std::unordered_map<int, int> b{
		{1, -1},
		{2, -2},
		{3, -3}
	};
	print(a);
	print(b);

	auto entries = split_string(" a");
	print(entries);

	ReactionNetwork rn("./decays_PDG21_massorder.dat");
	return 0;
}