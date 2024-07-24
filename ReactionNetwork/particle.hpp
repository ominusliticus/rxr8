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
struct Particle {

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
	void update(double dt, double temperature, RK4Stage stage);

	/// @brief Combine RK4 stages to complete single time step
	void finish_time_step();

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

inline void
Particle::finish_time_step()
{
	density += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6, 0;
	k1 = k2 = k3 = k4 = 0.0;
	already_visited   = false;
}

inline void
Particle::update(double dt, double temperature, RK4Stage stage)
{
	// if (!b_eq_density_set)
	// {
	// 	eq_density = gauss_quad(
	// 	    [&](double q) -> double
	// 	    {
	// 		    double energy{ 0.0 };
	// 		    if (std::fabs(temperature / mass) < 1e-2) energy = q * q / (2.0 * mass);
	// 		    else energy = std::sqrt(q * q + mass * mass);

	// 		    double f{ 0.0 };
	// 		    switch (spin_type)
	// 		    {
	// 			    case ParticleSpin::BOLTZ :
	// 			    {
	// 				    double thermal_wavelength{ std::sqrt(2.0 * pi / (mass * temperature)) };
	// 				    double norm{ 1.0 };
	// 				    for (auto i{ 0 }; i < 3; ++i)
	// 					    norm *= thermal_wavelength;
	// 				    f = std::exp(energy / temperature) / norm;
	// 				    break;
	// 			    }
	// 			    case ParticleSpin::FERMI :
	// 			    {
	// 				    f = 1.0 / (std::exp(energy / temperature) + 1.0);
	// 				    break;
	// 			    }
	// 			    case ParticleSpin::BOSON :
	// 			    {
	// 				    f = 1.0 / (std::exp(energy / temperature) - 1.0);
	// 				    break;
	// 			    }
	// 		    }

	// 		    return degeneracy * f / (2.0 * pi * pi) / (hbar * hbar * hbar);
	// 	    },
	// 	    0.0,
	// 	    inf,
	// 	    1e-10,
	// 	    3
	// 	);
	// 	b_eq_density_set = true;
	// }
	// auto decays = decay_width * density;
	// /// TODO: Need to continue editting here
	// auto inv_decays = 0.0;
	// update(-decays, dt, stage);
	// for (auto& info : reaction_infos)
	// 	info.propagate(decays, dt, stage);
	double *k;
	switch (stage)
	{
		case RK4Stage::FIRST :
			k = &k1;
			break;
		case RK4Stage::SECOND :
			k = &k2;
			break;
		case RK4Stage::THIRD :
			k = &k3;
			break;
		case RK4Stage::FOURTH :
			k = &k4;
			break;
	}
}
