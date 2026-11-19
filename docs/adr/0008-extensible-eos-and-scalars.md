# Extensible equation of state and scalar fields

Date: 2026-11-19. Status: accepted.

An equation of state is an immutable `isbits` struct with four methods: `pressure(eos, ρ)`, `enthalpy(eos, ρ)`, `sound_speed²(eos, ρ)` and `density_from_enthalpy(eos, h)`. `Isothermal` and `Polytropic` are built in. A user adds an EOS with one struct and four one-line methods, and these run inside GPU kernels without change.

The state stores its fields in a `NamedTuple` (`(:ux, :uy)` or `(:ux, :uy, :rho)`), so a future coupled scalar is one more named field plus its tendency. Extra scalars are not part of the current scope.
