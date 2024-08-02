#include "particle.hpp"

Particle::Particle(
    long        pid,
    double      mass,
    double      degeneracy,
    double      decay_width,
    SpinStat    spin_stat,
    std::size_t decay_channels
)
{
	m_pid         = pid;
	m_mass        = mass;
	m_degeneracy  = degeneracy;
	m_decay_width = decay_width;
	m_spin_stat   = spin_stat;
	m_reaction_infos.reserve(decay_channels);
}

void
Particle::finalize_time_step(void)
{
	m_density += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
	k1 = k2 = k3 = k4 = 0.0;
	m_already_visited = false;
}

void
Particle::update(double delta_density, double dt, RK4Stage stage)
{
	switch (stage)
	{
		case RK4Stage::FIRST :
		{
			k1 += dt * delta_density;
			break;
		}
		case RK4Stage::SECOND :
		{
			k2 += 0.5 * dt * delta_density;
			break;
		}
		case RK4Stage::THIRD :
		{
			k3 += 0.5 * dt * delta_density;
			break;
		}
		case RK4Stage::FOURTH :
		{
			k4 += dt * delta_density;
			break;
		}
	}
}

double
Particle::get_eq_density(double temperature)
{
	if (m_already_visited) return m_eq_density;

	m_eq_density = gauss_quad(
	    [&](double q) -> double
	    {
		    double energy{ 0.0 };
		    if (std::fabs(temperature / m_mass) < 1e-2) energy = q * q / (2.0 * m_mass);
		    else energy = std::sqrt(q * q + m_mass * m_mass);

		    double f{ 0.0 };
		    switch (m_spin_stat)
		    {
			    case SpinStat::MB :
			    {
				    double thermal_wavelength{ std::sqrt(2.0 * pi / (m_mass * temperature)) };
				    double norm{ 1.0 };
				    for (auto i{ 0 }; i < 3; ++i)
					    norm *= thermal_wavelength;
				    f = std::exp(energy / temperature) / norm;
				    break;
			    }
			    case SpinStat::FD :
			    {
				    f = 1.0 / (std::exp(energy / temperature) + 1.0);
				    break;
			    }
			    case SpinStat::BE :
			    {
				    f = 1.0 / (std::exp(energy / temperature) - 1.0);
				    break;
			    }
		    }

		    // Return density in units fm^{-3}
		    return m_degeneracy * f / (2.0 * pi * pi) / (hbar * hbar * hbar);
	    },
	    0.0,
	    inf,
	    1e-10,
	    3
	);
	return 0.0;
}

void
Particle::add_reaction(ReactionInfo&& info)
{
	m_reaction_infos.push_back(std::move(info));
}