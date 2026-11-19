#!/usr/bin/env bash
# r1-E7: a Julia package extension loads only when ALL its trigger packages are loaded.
# Package A declares AExt = ["B", "C"]. Loading B alone does not load AExt.
# So one NS2DMakieExt cannot trigger on "GLMakie or CairoMakie" (ADR-0012, ROADMAP:36,55).
# Run: bash r1-E7.sh   (offline; makes throwaway packages in a temp dir)
set -e
D=$(mktemp -d); cd "$D"
uid() { python3 -c 'import uuid;print(uuid.uuid4())'; }
for P in B C; do mkdir -p $P/src; printf 'name = "%s"\nuuid = "%s"\nversion = "0.1.0"\n' $P $(uid) > $P/Project.toml; echo "module $P end" > $P/src/$P.jl; done
BU=$(grep uuid B/Project.toml|cut -d'"' -f2); CU=$(grep uuid C/Project.toml|cut -d'"' -f2)
mkdir -p A/src A/ext
printf 'name = "A"\nuuid = "%s"\nversion = "0.1.0"\n[weakdeps]\nB = "%s"\nC = "%s"\n[extensions]\nAExt = ["B", "C"]\n' $(uid) $BU $CU > A/Project.toml
echo 'module A end' > A/src/A.jl; echo 'module AExt; using A, B, C; end' > A/ext/AExt.jl
JULIA_PKG_OFFLINE=true julia --startup-file=no --project=env -e 'using Pkg; Pkg.develop([PackageSpec(path="A"),PackageSpec(path="B"),PackageSpec(path="C")])' >/dev/null 2>&1
julia --startup-file=no --project=env -e 'using A, B; println("only B loaded: AExt loaded = ", Base.get_extension(A, :AExt) !== nothing); using C; println("B and C loaded: AExt loaded = ", Base.get_extension(A, :AExt) !== nothing)' 2>/dev/null
# Output on julia 1.12.4: "only B loaded: AExt loaded = false", "B and C loaded: AExt loaded = true"
