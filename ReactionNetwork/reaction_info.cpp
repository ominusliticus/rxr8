#include "reaction_info.hpp"
#include "particle.hpp"

void
ReactionInfo::calculate(std::shared_ptr<Particle> particle, double dt, double temperature, RK4Stage stage)
{
	switch (reaction_type)
	{
		case ReactionType::DECAY :
		{
			// Contribution from decays
			auto eq_density = particle->get_eq_density(temperature);
			auto density    = particle->get_density();
			density += particle->get_RK4Stage_offset(stage);
			auto from_decays = density / eq_density;

			// Contribution form inverse decays
			auto from_inv_decays = 1.0;
			for (auto product : products)
				from_inv_decays *= product->get_density() / product->get_eq_density(temperature);

			// Combine using dn/dt = Gamma n_eq (-n / n_eq + n1 n2 / (n1_eq n2_eq))
			// Which updates the current particles abundance
			auto delta_density = reaction_rate * eq_density * (from_inv_decays - from_decays);

			// update rk4 stage
			particle->update(delta_density, dt, stage);
			for (auto product : products)
				product->update(-delta_density, dt, stage);
		}
	}
}

void
ReactionInfo::calculate(std::shared_ptr<Particle> particle, double dt, double temperature, RK4Stage stage) const
{
	switch (reaction_type)
	{
		case ReactionType::DECAY :
		{
			// Contribution from decays
			auto eq_density = particle->get_eq_density(temperature);
			auto density    = particle->get_density();
			density += particle->get_RK4Stage_offset(stage);
			auto from_decays = density / eq_density;

			// Contribution form inverse decays
			auto from_inv_decays = 1.0;
			for (auto product : products)
				from_inv_decays *= product->get_density() / product->get_eq_density(temperature);

			// Combine using dn/dt = Gamma n_eq (-n / n_eq + n1 n2 / (n1_eq n2_eq))
			// Which updates the current particles abundance
			auto delta_density = reaction_rate * eq_density * (from_inv_decays - from_decays);

			// update rk4 stage
			particle->update(delta_density, dt, stage);
			for (auto product : products)
				product->update(-delta_density, dt, stage);
		}
	}
}