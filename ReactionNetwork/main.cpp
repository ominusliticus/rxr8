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

	std::string_view data_sheet{ "../input/PDG21Plus/hadron_lists/PDG21Plus/PDG21Plus_massorder.dat" };
	std::string_view decay_sheet{ "../input/PDG21Plus/hadron_lists/PDG21Plus/full_decays/decays_PDG21_massorder.dat" };
	ReactionNetwork  rn(data_sheet, decay_sheet);
	return 0;
}