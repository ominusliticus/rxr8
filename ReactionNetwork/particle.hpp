#pragma once

#include <cstddef>
#include <cstdio>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "../constants.hpp"
#include "../integration.hpp"

#include "print.hpp"
#include "reaction_info.hpp"
#include "reaction_type.hpp"
#include "rk4_stages.hpp"
#include "spin_statistics.hpp"

/// @brief Stores the particle ID and reactions, and facilitates density updates
/// @details Stores the particle ID, (for now, only the) decay width, and list of daughters with the corresponding
/// branching ratio (will contain other reaction rates in the future). Daughters are stored within the `ReactionInfo`
/// structure, and just consists of a list of shared pointers to the daughter particles. Here the design choice was
/// such that the `Node` class did not have to be aware of the particle dictionary in the `ReactionNetwork` class.
/// The class also contains facilities that interoperate with the `ReactionNetwork` class to perform time-stepping
/// in a self-consistent way: that is, the first stage of of Runge-Kutta should be completed before the second starts
/// and so on.
class Particle
{
	public:
	Particle() = default;

	Particle(
	    long        pid,
	    double      mass,
	    double      degeneracy,
	    double      decay_width,
	    SpinStat    spin_stat,
	    std::size_t decay_channels
	);

	double get_density(void) { return m_density; }

	void set_density(double density) { m_density = density; }

	int get_pid(void) { return m_pid; }

	void   update(double delta_density, double dt, RK4Stage stage);
	void   finalize_time_step(void);
	double get_eq_density(double temperature);
	double get_RK4Stage_offset(RK4Stage stage);
	void   add_reaction(ReactionInfo&& info);

	std::vector<ReactionInfo> const& get_reactions(void) const { return m_reaction_infos; }

	public:
	// Stages for 4th-order Runge-Kutta
	double k1{ 0.0 };    // First stage of RK4
	double k2{ 0.0 };    // Second stage of RK4
	double k3{ 0.0 };    // Third stage of RK4
	double k4{ 0.0 };    // Fourth stage of RK4

	private:
	// For reactions
	SpinStat                  m_spin_stat;
	long                      m_pid;
	double                    m_eq_density;
	double                    m_density;
	double                    m_mass;
	double                    m_decay_width;
	double                    m_degeneracy;
	std::vector<ReactionInfo> m_reaction_infos;
	bool                      m_eq_density_calculated{ false };
};