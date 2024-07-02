#pragma once

#include <cassert>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>

#include "print.hpp"
#include "string_utility.hpp"

#include "node.hpp"

// Ideas:
//  - Looping over span tree can be made parallel by having multiple threads loop over tree and making input variable an
//  atomic.
//  - Needs an detailed balance contribution for each reaction (i.e., every decay should have its corresponding inverse
//  dacay)

/// @brief Structure that stores and evolves the densities of particles
/// @details This class provides the functionality that stores a list of particles, their initial densities and then
/// integrates their rate equations in time using a Runge-Kutta 4th order time-stepping scheme.
class ReactionNetwork
{
	public:
	ReactionNetwork() = default;
	ReactionNetwork(std::string_view decays_file);

	/// TODO: Add version that takes a dictionary and another that has an interable container
	void set_initial(std::string const& file);    // set all initial values
	void time_step(double dt, RK4Stage stage);    // Perform a single time step, involves traversing spanning tree
	void finish_time_step();                      // Update densities based on input and set input to zero

	private:
	std::unordered_map<long, std::shared_ptr<Node>> m_dict;
};

/// @brief Constructor for structure that stores and evolves the densities of particles
/// @param decays_file path to file storing the reaction information sheet in mass-ordering form
/// @details This class provides the functionality that stores a list of particles, their initial densities and then
/// integrates their rate equations in time using a Runge-Kutta 4th order time-stepping scheme. Function can fail due to
/// file not exsint, and will terminate program
inline ReactionNetwork::ReactionNetwork(std::string_view decays_file)
{
	std::fstream fin(decays_file.data(), std::fstream::in);
	assert(fin.is_open() && "Decays file failed to open");
	std::string line;
	// File layout (by column name) [all units in GeV]
	// PID Name Mass Width Spin B S Q C B I Iz Q No.-decays
	// PID No.-daughters Branching-ratio PID-1 PID-2 PID-2 PID-4 PID-5

	// Read in massorder file and parser into reaction network
	long first_pid{ 0 };
	// int  counter{ 0 };
	while (!fin.eof())
	{
		std::getline(fin, line);
		auto entries{ split_string(line) };
		auto pid{ std::stol(entries[0]) };
		auto width{ std::stod(entries[3]) };
		auto num_decays{ std::stoi(entries.back()) };

		if (!first_pid) first_pid = pid;

		if (m_dict.contains(pid)) m_dict[pid]->decay_width = width;
		else m_dict[pid] = std::make_shared<Node>(pid, width);
		for (int i{ 0 }; i < num_decays; ++i)
		{
			std::getline(fin, line);
			auto entries{ split_string(line) };
			auto n_daughters{ std::stoi(entries[1]) };
			auto br{ std::stod(entries[2]) };

			std::vector<std::shared_ptr<Node>> daughters;
			for (int n = 0; n < n_daughters; ++n)
			{
				auto daughter_pid{ std::stol(entries[4 + n]) };
				if (m_dict.contains(daughter_pid)) daughters.push_back(m_dict[daughter_pid]);
				else
				{
					m_dict[daughter_pid] = std::make_shared<Node>(daughter_pid, 0.0);
					daughters.push_back(m_dict[daughter_pid]);
				}
			}
			Node::ReactionInfo ri{ .branching_ratio = br, .decay_products = std::move(daughters) };
			m_dict[pid]->reaction_infos.push_back(std::move(ri));
		}
	}
	// build_minimum_spanning_tree(m_dict[first_pid]);
}

/// @brief Takes an initial condition file and initializes all the density for the reaction newtork
/// @details Function can fail, and results in termination of the program
/// @param file std::string of file containing initial conditon file
inline void
ReactionNetwork::set_initial(std::string const& file)
{
	std::fstream fin(file.data(), std::fstream::in);
	assert(fin.is_open() && "Initial condition file failed to open");
	long   pid;
	double density;

	while (!fin.eof())
	{
		fin >> pid >> density;
		m_dict[pid]->density = density;
	}
}

/// @brief Preforms a partial time integration step of the Runge-Kutta 4th order algorithm
/// @param double dt size of single time time step
/// @param RK4Stage value from the enum class indicating which stage in the Runge-Kutta fourth order scheme to perform
inline void
ReactionNetwork::time_step(double dt, RK4Stage stage)
{
	for (auto [key, particle] : m_dict)
		particle->propagate(dt, stage);
}

/// @brief Combine the inidividual Runge-Kutte 4th order stages to preform update of particle densities after one full
/// time step
inline void
ReactionNetwork::finish_time_step()
{
	for (auto [key, particle] : m_dict)
		particle->finish_time_step();
}