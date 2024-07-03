
#include <memory>
#include <string>
#include <vector>

/// @brief Enum class that conveniently indicates which Runge-Kutta 4th order stage to perform
enum class RK4Stage { FIRST, SECOND, THIRD, FOURTH };

/// @brief Stores the particle ID and reactions, and facilitates density updates
/// @details Stores the particle ID, (for now, only the) decay width, and list of daughters with the corresponding
/// branching ratio (will contain other reaction rates in the future). Daughters are stored within the `ReactionInfo`
/// structure, and just consists of a list of shared pointers to the daughter particles. Here the design choice was
/// such that the `Node` class did not have to be aware of the particle dictionary in the `ReactionNetwork` class.
/// The class also contains facilities that interoperate with the `ReactionNetwork` class to perform time-stepping
/// in a self-consistent way: that is, the first stage of of Runge-Kutta should be completed before the second starts
/// and so on.
struct Node {

	Node() = default;

	Node(long pid_, double width)
	{
		pid         = pid_;
		decay_width = width;
	}

	void update(double delta_density, double dt, RK4Stage stage);
	void finish_time_step();
	void propagate(double dt, RK4Stage stage);

	struct ReactionInfo {
		enum class Type { DECAY };
		void propagate(double decays, double dt, RK4Stage stage);

		double                             branching_ratio;
		std::vector<std::shared_ptr<Node>> decay_products;
	};

	// For reactions
	long                      pid;
	double                    density;
	double                    decay_width;
	double                    k1{ 0.0 };    // First stage of RK4
	double                    k2{ 0.0 };    // Second stage of RK4
	double                    k3{ 0.0 };    // Third stage of RK4
	double                    k4{ 0.0 };    // Fourth stage of RK4
	std::vector<ReactionInfo> reaction_infos;
	// For minimum spanning tree
	std::vector<std::shared_ptr<Node>> spanning_nodes;
	bool                               already_visited{ false };
};

// Utility functions
inline void
Node::update(double delta_density, double dt, RK4Stage stage)
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
			k2 = dt * (delta_density + k1 / 2.0);
			break;
		}
		case RK4Stage::THIRD :
		{
			k3 = dt * (delta_density + k2 / 2.0);
			break;
		}
		case RK4Stage::FOURTH :
		{
			k4 = dt * (delta_density + k3);
			break;
		}
	}
}

inline void
Node::finish_time_step()
{
	density += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6, 0;
	k1 = k2 = k3 = k4 = 0.0;
}

inline void
Node::propagate(double dt, RK4Stage stage)
{
	auto decays = decay_width * density;
	update(-decays, dt, stage);
	for (auto& info : reaction_infos)
		info.propagate(decays, dt, stage);
}

inline void
Node::ReactionInfo::propagate(double decays, double dt, RK4Stage stage)
{
	for (auto node : decay_products)
		node->update(branching_ratio * decays, dt, stage);
}
