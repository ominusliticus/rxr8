#include "particle.hpp"

inline void
Particle::finish_time_step(void)
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