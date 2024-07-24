# `ReactionType` enumeration class

Contains the type of reactions considered in the the reaction network.
Reactiont types that can be considered, but have not been included yet `TWO_TO_TWO`, `RESONANCE`, `THREE_TO_TWO` or `TO_TO_THREE` are possible extensions to this `enum class`.
This enumeration is used to customize how the update to the densities are calculated.

## Entries

- `DECAY`

<!-- ==================================================================== -->

# `RK4Stage` enumeration class

Contains enumerations for the four stages of fourth-order Runge-Kutta time integration scheme.
Used while time-stepping through reaction network.

## Entries
- `FIRST`
- `SECOND`
- `THIRD`
- `FOURTH`

<!-- ==================================================================== -->

# `SpinStat` enumeration class

Used to determine what particle distribution needs to be used to calculate the equilibrium densities

## Entries

- `MB`: Maxwell-Boltzmann distribution
- `FD`: Fermi-Dirac distribution
- `BE`: Bose-Einstein distribution

<!-- ==================================================================== -->

# `ReactionInfo` class

This class stores reaction parameters, reactants, and products for a given reaction that a particle can undergo.
Every particle, owns a list of `ReactionInfo`s that is looped over during time stepping.
This information is used to calculate the updated density for the particle owning the instance of `ReactionInfo` and all its products in the reaction.

## Member variables

- `reaction_type`: (`ReactionType`) allows for customizable behavior on how to calculate the reaction rate. 
- `reaction_rate`: (`double`) the parameters that controls how much of the density is gained/lost at each time step. (For decays, this is the decay width times the branching fraction)
- `reactants`: (`std::vector<std::shared_ptr<Particle>>`) vector of other particles participating in the reaction (leads to some duplicated calculations, optimizations should be considered)
- `products`: (`std::vector<std::shared_ptr<Particle>>`) vector of products of the reaction, which have their densities updated 

## Member functions

### `ReactionInfo::calculate`

Calculates the `delta_density` given the time step `dt` and background temperature `temperature`, then subtracts `delta_density` from the particle _owning_ this instance of `ReactionInfo` and adding `delta_density` to the particles in `products`.

#### Signature and return value
```c++
calculate(double density, double eq_density, double dt, double temperature) -> void
```

#### Function parameters

- `density`: (`double`) the density of the particle _owning_ this instance of the `ReactionInfo` struct
- `eq_density`: (`double`) the equilibrium density, given temperature `temperature`, for the particle _owning_ this instance of `ReactionInfo`
- `dt`: (`double`) the time step size
- `temperature`: (`double`) the background temperature

<!-- ==================================================================== -->

# `Particle` class

The main ingredient in the reaction network.
Each instance represents a unique particle (pid), and stores the information about all its possible reactants and products, and the reaction strength

## Member variables

- `pid`: (`long`) the unique identifier assigned to each particle by particle physicists
- `spin_stat`: (`SpinStat`) determines what distribution to use to calculate equilibrium density
- `mass`: (`double`) mass of the particle; needed to calculate the equilibrium density
- `temperature`: (`double`) the background temperature used to calculate the equilibrium density
- `degeneracy`: (`double`) the spin, isospin, and other internal d.o.f. degeneracy for the particle being considered
- `density`: (`double`) the density of the particle being evolved
- `eq_density`: (`double`) the member variable that stores the equilibrium density for every time step
- `already_visited`: (`bool`) a variable that is reset at the end of every time step and keeps track of which equilibrium densities have already been calculated while iterating through the particle list
- `reaction_info`: (`std::vector<ReactionInfo>`) list of `ReactionInfo` instances that are used to calculate the density updates
- `k1`,`k2`,`k3`,`k4`: (`double`) {initialized to zero} stores the update values from each stage of the fourth-order Runge-Kutta scheme

## Member functions

### `Particle` constructors

Takes information from particle data sheets and converts into an particle instance

#### Signature and return value

```c++
Particle(long pid_, double mass_, double degeneracy_, SpinStat spin_stat_)
```

#### Function parameters

### `Particle::update`

Receives the `delta_density` calculated from whatever reaction was considered and adds it to respective variable (`k1`, `k2`, `k3`, or `k4`) corresponding to stage `stage`

#### Signature and return value

```c++
update(double delta_density, double dt, RK4Stage stage) -> void
```
#### Function parameters

### `Particle::finalize_time_step`

Combines the four stages and sets `k1` through `k4` to zero for the next time step.
Also marks `already_visited` as false so the equilibrium densities can be calculated on the next time step.

#### Signature and return value

```c++
finalize_time_step(void) -> void
```

#### Function parameters