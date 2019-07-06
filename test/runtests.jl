using PyPlot
using JLD2, FileIO
using PumpProbeFitting
using PumpProbeModels

data = load("/Users/lackner/Documents/DataSFG/Notebooks/CaAra/export/pp_d+fr_dif.jld2")

x = data["wn3"]
y = data["dl"]
z = data["mdiff3[roi...]"]

m = PumpProbeModel(x, y, 2, pumpmode=:δ)
m.parameters.A = [4.9, 2.1]
m.parameters.Γ = [7.0, 8.0]
m.parameters.ω = [2860.0, 2935.0]
m.parameters.τ[1] = [90.0, 200.0]
m.parameters.σ[1] = [0.0, 0.0]
m.parameters.γ[1] = [0.0, 0.0]

lower, upper = PumpProbeFitting.generate_bounds(m)

r,g = PumpProbeFitting.fit(z, m, lower, upper)

F = Array{Float64}(undef, length(z[:]))
xx = eltype(m.wavenumbers)[]
tt = eltype(m.delaytimes)[]
for _x in m.wavenumbers, _t in m.delaytimes
    push!(xx, _x)
    push!(tt, _t)
end
evaluate!(F, xx, tt, m)
FR = reshape(F, size(z))

figure(figsize=(9,4))
subplot(121)
pcolormesh(z)
subplot(122)
pcolormesh(FR)
