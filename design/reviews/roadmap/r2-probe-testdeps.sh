#!/usr/bin/env bash
# r2-probe-testdeps: a package lists FFTW only in [extras] + [targets] (as WP-0, ROADMAP:91, plans
# for CUDA and AMDGPU). The commands of CLAUDE.md:16-18 run test scripts with `julia --project`,
# outside Pkg.test. There an [extras] package does not load. JULIA_LOAD_PATH="@:@stdlib" removes
# the user's default environment, as on a machine where AMDGPU is not installed globally.
# Run: bash r2-probe-testdeps.sh   (offline; FFTW must be in the depot)
set -e
D=$(mktemp -d); mkdir -p "$D/P/src" "$D/P/test"
cat > "$D/P/Project.toml" <<TOML
name = "P"
uuid = "5b1e4c1e-0000-4000-8000-000000000001"
version = "0.1.0"
[extras]
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
[targets]
test = ["Test", "FFTW"]
TOML
echo 'module P end' > "$D/P/src/P.jl"
echo 'using Test, FFTW; println("test deps loaded")' > "$D/P/test/runtests.jl"
cd "$D/P"
JULIA_PKG_OFFLINE=true julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()' >/dev/null 2>&1 || true
echo "--- julia --project test/runtests.jl (CLAUDE.md:16 style):"
JULIA_LOAD_PATH="@:@stdlib" julia --startup-file=no --project=. test/runtests.jl 2>&1 | grep -m1 -iE "error|loaded" || true
