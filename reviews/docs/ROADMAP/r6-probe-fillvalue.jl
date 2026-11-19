# r6-probe-fillvalue: ROADMAP:146 and ADR-0018:9 call NaN "the NetCDF fill value". NCDatasets writes NaN as a value,
# sets no _FillValue, and reads it back as NaN (not missing); the NetCDF default fill for double is 9.97e36.
# Needs NCDatasets: julia -e 'using Pkg; Pkg.activate(temp=true); Pkg.add("NCDatasets"); include("r6-probe-fillvalue.jl")'
using NCDatasets
f = tempname() * ".nc"; isfile(f) && rm(f)
NCDataset(f, "c") do ds
    defDim(ds, "time", Inf)
    v = defVar(ds, "inj_E_kin", Float64, ("time",))
    v[1:2] = [NaN, 1.0]
end
NCDataset(f) do ds
    v = ds["inj_E_kin"]
    println("attribs: ", Dict(v.attrib))
    println("read (CF-decoded): ", v[:])
    println("NetCDF default fill for double: ", NCDatasets.fillvalue(Float64))
end
