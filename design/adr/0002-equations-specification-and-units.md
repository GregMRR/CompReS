# Equations specification and units

Date: 2026-11-19. Status: accepted.

`docs/ns2d_equations.tex` (v2, 2026-11-19) is the binding specification of the equations, invariants, spectra and fluxes. A change to an equation is a new version of that document plus a new ADR. The derivations of v2 were checked by hand during the design phase; no error was found.

Units: the mean density is `ρ₀ = ⟨ρ⟩ = 1` in both models. Mass is conserved, so this only fixes the density unit. With it, `E_kin` and `Z` have the same meaning in both models, and `Ω = Z` in the incompressible model.

Both models use the velocity formulation; the code never evolves vorticity or momentum density. The compressible state is `(u, ρ)`.
