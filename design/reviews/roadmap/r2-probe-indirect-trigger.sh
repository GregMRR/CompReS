#!/usr/bin/env bash
# r2-probe-indirect-trigger: does an extension with trigger B load when the user loads only C,
# a package that depends on B (as GLMakie depends on Makie)? Checks the ADR-0012 amendment.
# Run: bash r2-probe-indirect-trigger.sh   (offline; throwaway packages in a temp dir)
set -e
D=$(mktemp -d); cd "$D"
uid() { python3 -c 'import uuid;print(uuid.uuid4())'; }
mkdir -p B/src; BU=$(uid); printf 'name = "B"\nuuid = "%s"\nversion = "0.1.0"\n' $BU > B/Project.toml; echo "module B end" > B/src/B.jl
mkdir -p C/src; printf 'name = "C"\nuuid = "%s"\nversion = "0.1.0"\n[deps]\nB = "%s"\n' $(uid) $BU > C/Project.toml; echo "module C; using B; end" > C/src/C.jl
mkdir -p A/src A/ext
printf 'name = "A"\nuuid = "%s"\nversion = "0.1.0"\n[weakdeps]\nB = "%s"\n[extensions]\nAExt = "B"\n' $(uid) $BU > A/Project.toml
echo 'module A end' > A/src/A.jl; echo 'module AExt; using A, B; end' > A/ext/AExt.jl
JULIA_PKG_OFFLINE=true julia --startup-file=no --project=env -e 'using Pkg; Pkg.develop([PackageSpec(path="B"), PackageSpec(path="A"), PackageSpec(path="C")]); Pkg.rm("B")' >/dev/null 2>&1
echo "env direct deps: $(grep -A5 '\[deps\]' env/Project.toml | grep -oE '^[A-Z]' | tr '\n' ' ')"
julia --startup-file=no --project=env -e 'using A, C; println("only C loaded (B indirect): AExt loaded = ", Base.get_extension(A, :AExt) !== nothing)' 2>&1 | tail -1
