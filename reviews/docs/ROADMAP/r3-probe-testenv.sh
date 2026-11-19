#!/usr/bin/env bash
# r3-probe-testenv: ROADMAP:42, 92 and CLAUDE.md:16-18 use test/Project.toml with NS2D through
# [sources] path "..". Check on a toy package that Pkg.test() and `julia --project=test test/runtests.jl`
# both work, and that `--project=docs` has an environment (ROADMAP:67 lists no docs/Project.toml).
# Run: bash r3-probe-testenv.sh <scratch dir>
set -u; D=${1:?scratch dir}/testenv; rm -rf "$D"; mkdir -p "$D/src" "$D/test"
cat > "$D/Project.toml" <<'T'
name = "Toy"
uuid = "1d6a5b0e-0000-4000-8000-000000000001"
version = "0.1.0"
T
echo 'module Toy; f() = 1; end' > "$D/src/Toy.jl"
cat > "$D/test/Project.toml" <<'T'
[deps]
Toy = "1d6a5b0e-0000-4000-8000-000000000001"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
[sources]
Toy = {path = ".."}
T
echo 'using Toy, Test; @test Toy.f() == 1; println("runtests ok")' > "$D/test/runtests.jl"
cd "$D"
export JULIA_PKG_OFFLINE=true JULIA_DEPOT_PATH="$D/depot:"
echo "--- Pkg.test()";  julia --project -e 'using Pkg; Pkg.test()' 2>&1 | grep -E 'runtests ok|ERROR|tests passed' | head -5
echo "--- julia --project=test test/runtests.jl"; julia --project=test -e 'using Pkg; Pkg.instantiate()' >/dev/null 2>&1; julia --project=test test/runtests.jl 2>&1 | tail -2
