#pragma once

#include <cassert>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>

#include "print.hpp"
#include "reaction_type.hpp"
#include "rk4_stages.hpp"
#include "string_utility.hpp"

#include "particle.hpp"
#include "reaction_info.hpp"

/// @brief Structure that stores and evolves the densities of particles
/// @details This class provides the functionality that stores a list of particles, their initial densities and then
/// integrates their rate equations in time using a Runge-Kutta 4th order time-stepping scheme.
class ReactionNetwork
{
	public:
	ReactionNetwork() = default;
	ReactionNetwork(std::string_view particle_datasheet, std::string_view particle_reactions);

	void time_step(double dt, double temperature);
	void finalize_time_step();

	private:
	std::unordered_map<long, std::shared_ptr<Particle>> m_particles;
};

/// @brief Constructor for structure that stores and evolves the densities of particles
/// @param decays_file path to file storing the reaction information sheet in mass-ordering form
/// @details This class provides the functionality that stores a list of particles, their initial densities and then
/// integrates their rate equations in time using a Runge-Kutta 4th order time-stepping scheme. Function can fail due to
/// file not exist, and will terminate program
inline ReactionNetwork::ReactionNetwork(std::string_view particle_datasheet, std::string_view particle_deacys)
{
	std::fstream fin(particle_datasheet.data(), std::fstream::in);
	assert(fin.is_open() && "Reactions file failed to open");
	std::string line;
	// File layout (by column name) [all units in GeV]
	// PID Name Mass Width Spin-Degen. B S c b I Iz Q Num-decays
	while (!fin.eof())
	{
		std::getline(fin, line);
		auto entries{ split_string(line) };
		auto pid{ std::stol(entries[0]) };
		auto mass{ std::stod(entries[2]) };
		auto width{ std::stod(entries[3]) };
		auto spin_degen{ std::stod(entries[4]) };
		auto num_decays{ std::stoull(entries.back()) };
		auto spin_stat{ static_cast<int>(spin_degen) % 2 == 0 ? SpinStat::FD : SpinStat::BE };

		m_particles[pid] = std::make_shared<Particle>(pid, mass, spin_degen, width, spin_stat, num_decays);
	}
	fin.close();

	fin.open(particle_deacys.data(), std::fstream::in);
	assert(fin.is_open() && "Reactions file failed to open");
	// File layout (by column name) [all units in GeV]
	// PID Name Mass Width Spin B S Q C B I Iz Q No.-decays
	// PID No.-daughters Branching-ratio PID-1 PID-2 PID-2 PID-4 PID-5

	while (!fin.eof())
	{
		std::getline(fin, line);
		auto entries{ split_string(line) };
		auto pid{ std::stol(entries[0]) };
		auto width{ std::stod(entries[3]) };
		auto num_decays{ std::stoi(entries.back()) };

		for (int i{ 0 }; i < num_decays; ++i)
		{
			std::getline(fin, line);
			auto entries{ split_string(line) };
			auto n_daughters{ std::stoi(entries[1]) };
			auto br{ std::stod(entries[2]) };

			std::vector<std::shared_ptr<Particle>> reactants{ m_particles[pid] };
			std::vector<std::shared_ptr<Particle>> products;
			for (int n = 0; n < n_daughters; ++n)
				products.push_back(m_particles[std::stol(entries[4 + n])]);
			ReactionInfo ri{ .reaction_type = ReactionType::DECAY,
				             .reaction_rate = br * width,
				             .reactants     = std::move(reactants),
				             .products      = std::move(products) };
			m_particles[pid]->add_reaction(std::move(ri));
		}
	}
	// build_minimum_spanning_tree(m_dict[first_pid]);
}

/// @brief Preforms a partial time integration step of the Runge-Kutta 4th order algorithm
/// @param double dt size of single time time step
/// @param RK4Stage value from the enum class indicating which stage in the Runge-Kutta fourth order scheme to perform
inline void
ReactionNetwork::time_step(double dt, double temperature)
{
	for (auto stage : std::vector<RK4Stage>{ RK4Stage::FIRST, RK4Stage::SECOND, RK4Stage::THIRD, RK4Stage::FOURTH })
		for (auto [key, particle] : m_particles)
			for (auto const& reaction : particle->get_reactions())
				reaction
				    .calculate(particle->get_density(), particle->get_eq_density(temperature), dt, temperature, stage);
}

/// @brief Combine the individual Runge-Kutte 4th order stages to preform update of particle densities after one full
/// time step
inline void
ReactionNetwork::finalize_time_step()
{
	for (auto [key, particle] : m_particles)
		particle->finalize_time_step();
}