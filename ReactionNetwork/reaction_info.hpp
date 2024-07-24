#pragma once

#include <memory>
#include <unistd.h>
#include <vector>

#include "reaction_type.hpp"
#include "rk4_stages.hpp"

class Particle;

/// @brief An internal struct to class the stores the reaction details for each class
/// @details This struct stores the reaction details for all process that have been supplied for a given particle
/// It also does the job of calculating the loss of the current particle for a given  time step, and propagating
/// this loss to the daughter particles of all the reactions
struct ReactionInfo {
	/// TODO: Needs a constructor

	/// @brief Calculates the amount of densities to pass along to all daughter particles
	/// @details Given a reaction `Type`, this function propogates the respective amount of this
	/// particle's density to all the products.
	/// @param dt double the current time step size
	/// @param temperature double the background temperature corresponding to the current time step
	/// @param stage RK4Stage enum class specifying which 4th order Runge-Kutta stage we are one
	/// @param type Type enum  class specifying what kind of reaction is to be considered
	/// @param args (variadic parameter) stores the necessary parameters to complete the calculation
	/// for a reaction of type `type`. The expected order of the parameters is as follows
	/// 	- Type::DECAY:
	///         Expected arguments: (current_density,)
	void calculate(double density, double eq_density, double dt, double temperature);

	ReactionType                           reaction_type;
	double                                 reaction_rate;
	std::vector<std::shared_ptr<Particle>> reactants;
	std::vector<std::shared_ptr<Particle>> products;
};

inline void
ReactionInfo::calculate(double density, double eq_density, double dt, double temperature)
{
}