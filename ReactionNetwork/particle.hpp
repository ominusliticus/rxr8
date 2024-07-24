#pragma once

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

	Particle(long pid_, double mass_, double degeneracy_, SpinStat spin_type_)
	{
		pid        = pid_;
		mass       = mass_;
		degeneracy = degeneracy_;
		spin_type  = spin_type_;
	}

	/// @brief Calculates the amount of densities to pass along to all daughter particles
	/// @details Given a reaction `Type`, this function propogates the respective amount of this
	/// particle's density to all the products.
	/// @param dt double the current time step size
	/// @param temperature double the background temperature corresponding to the current time step
	/// @param stage RK4Stage enum class specifying which 4th order Runge-Kutta stage we are one
	/// @param type ReactionType enum  class specifying what kind of reaction is to be considered
	/// @param args (variadic parameter) stores the necessary parameters to complete the calculation
	/// for a reaction of type `type`. The expected order of the parameters is as follows
	/// 	- Type::DECAY:
	///         Expected arguments: (current_density,)
	void update(double delta_density, double dt, RK4Stage stage);

	/// @brief Combine RK4 stages to complete single time step
	void finish_time_step(void);

	double get_eq_density(double temperature);
	void   add_reaction(ReactionInfo&& info);

	private:
	// For reactions
	SpinStat                  spin_type;
	long                      pid;
	double                    eq_density;
	double                    density;
	double                    mass;
	double                    degeneracy;
	double                    k1{ 0.0 };    // First stage of RK4
	double                    k2{ 0.0 };    // Second stage of RK4
	double                    k3{ 0.0 };    // Third stage of RK4
	double                    k4{ 0.0 };    // Fourth stage of RK4
	std::vector<ReactionInfo> reaction_infos;
	bool                      already_visited{ false };
};