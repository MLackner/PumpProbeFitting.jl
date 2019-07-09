using JLD2, FileIO
using PumpProbeFitting
using PumpProbeModels
using PumpProbePlotting
using PyPlot
using Images

data = load("/Users/lackner/Documents/DataSFG/Notebooks/CaAra/export/pp_d+fr_dif.jld2")
x = data["wn3"]
y = data["dl3"]
z = data["z3"]

xf = 0.2
yf = 0.5

NX = trunc(Int, length(x) * xf)
NY = trunc(Int, length(y) * yf)

x = imresize(x, NX)
y = imresize(y, NY)
z = imresize(z, NY, NX)

## START PARAMETERS
m = PumpProbeModel(x, y, 2, pumpmode=:triangle)
m.pumptimes = -19.0:4:9.0
m.parameters.A = [3.0, 3.0]
m.parameters.a = [0.1, 0.1]
m.parameters.ω = [2861.524, 2928.241]
m.parameters.σ[1] = [15.0, 15.0]
m.parameters.Δω = [0.0, 20.0]
m.parameters.t0 = 8.0
m.parameters.Γ = [9.0, 6.0]
m.parameters.τ[1] = [110.0, 160.0]
m.parameters.γ[1] = [0.0, 0.0]

## BOUNDS
lower, upper = PumpProbeFitting.generate_bounds(m)
lower.ω = [2860.0, 2927.0]
upper.ω = [2862.0, 2929.0]
lower.Δω = [-20.0, -5.0]
upper.Δω = [20.0, 30.0]
lower.σ[1] = [5.0, 5.0]
upper.σ[1] = [60.0, 60.0]
lower.t0 = 5.0
upper.t0 = 15.0
lower.τ[1] = [50.0, 120.0]
upper.τ[1] = [120.0, 250.0]
lower.A = [2.0, 2.0]
upper.A = [5.0, 5.0]
lower.a = [0.0, 0.0]
upper.a = [0.5, 0.5]
lower.Γ = [6.0, 5.0]
upper.Γ = [10.0, 8.0]

println("Fitting...")
@time r = fit(z, m_fit, lower, upper, maxIter=100)

printparams(m.parameters, lower, upper, r.m.parameters)

if r.f.converged == true
    @info "Fit Converged"
else
    @info "Fit Did Not Converge!"
end

println("RSQUARED: $(rsquared(r))")

dls = trunc.(Int, [20:5:60...] .* yf)
wns = trunc.(Int, [45, 65, 120, 140, 160] .* xf)

plot(r, m, dls, wns)
