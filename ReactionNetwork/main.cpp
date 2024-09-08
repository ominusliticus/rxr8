#include "particle.hpp"
#include "print.hpp"
#include "reaction_info.hpp"
#include "reaction_network.hpp"
#include "string_utility.hpp"

#include <filesystem>

double
ideal_hydro_temp(double tau, double tau_0, double T_0)
{
	return T_0 * std::exp(4.0 / 3.0 * std::log(tau_0 / tau));
}

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

	auto cwd{ std::filesystem::current_path() };
	auto hadron_list{ cwd / "../input/PDG21Plus/hadron_lists/PDG21Plus/PDG21Plus_massorder.dat" };
	auto decays_list{ cwd / "../input/PDG21Plus/hadron_lists/PDG21Plus/full_decays/decays_PDG21Plus_massorder.dat" };
	// auto             hadron_list{ cwd / "../test/particles.dat" };
	// auto             decays_list{ cwd / "../test/decays.dat" };
	std::string_view data_sheet{ hadron_list.c_str() };
	std::string_view decay_sheet{ decays_list.c_str() };
	ReactionNetwork  rn(data_sheet, decay_sheet);
	print(rn.get_particle_list()[111]->get_reactions().size());
	for (auto const& reaction : rn.get_particle_list()[111]->get_reactions())
		print("   ", reaction.products.size(), reaction.products[0]->get_pid(), reaction.reactants[0]->get_pid());

	double tau_0{ 0.1 };
	double dtau{ tau_0 / 20.0 };
	double tau_f{ 20.0 };
	double temperature{ 0.500 };

	rn.initialize_system(tau_0, temperature);
	for (auto tau = tau_0; tau <= tau_f; tau += dtau)
	{
		rn.time_step(dtau, ideal_hydro_temp(tau, tau_0, temperature));
		print(tau, rn.get_particle_density(111));
	}
	return 0;
}