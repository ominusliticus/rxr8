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

	double get_particle_density(long pid) { return m_particles[pid]->get_density(); }

	auto& get_particle_list() { return m_particles; }

	private:
	std::unordered_map<long, std::shared_ptr<Particle>> m_particles;
};