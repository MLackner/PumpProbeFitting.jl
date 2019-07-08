using JLD2, FileIO
using PumpProbeFitting
using PumpProbeModels
using Makie, AbstractPlotting

function plot(m::PumpProbeModel)
    x = m.wavenumbers
    y = m.delaytimes
    z = zeros(length(y), length(x))
    X, Y, Z = PumpProbeModels.convert1D(x, y, z)

    evaluate!(Z, X, Y, m)
    z .= reshape(Z, (length(y), length(x))) * -1

    scene = wireframe(y, x, z, color=(:red, 1.0))
    axis = scene[Axis]
    axis[:names, :axisnames] = ("Pump Delay (ps)", "Wavenumber (1/cm)", "")
    axis[:names, :rotation] = (1π,  π/2, π/2)
    axis[:names, :align] = ((:center, :center), (:center, :center), (:center, :center))
    axis[:names, :gap] = (10, 10)
    # scale!(scene, 1, 1, 500)
    # update_cam!(scene, FRect3D(Vec3f0(0), Vec3f0(1)))

    return scene
end

function plot(m::PumpProbeModel, data)
    scene = plot(m)
    z = deepcopy(data) .* -1
    scatter!(scene, m.delaytimes, m.wavenumbers, z, markersize=2, color=(:black, 0.5))

    # scale!(scene, 1, 1, 500)
    # update_cam!(scene, FRect3D(Vec3f0(0), Vec3f0(1)))

    scene
end

# data = load("/Users/lackner/Documents/DataSFG/Notebooks/CaAra/export/pp_d+fr_dif.jld2")

# x = data["wn3"]
# y = data["dl"]
# z = data["mdiff3[roi...]"]

x = range(2800.0, 3000.0, length=51)
y = range(-20, 300, length=61)
z = zeros(length(y), length(x))

X, Y, Z = PumpProbeModels.convert1D(x, y, z)

m = PumpProbeModel(x, y, 2, pumpmode=:triangle)
m.pumpfunction = y -> exp(-(y)^2/200)
m.parameters.A = [4.9, 5.1]
m.parameters.Γ = [7.0, 8.0]
m.parameters.ω = [2860.0, 2935.0]
m.parameters.τ[1] = [90.0, 200.0]
m.parameters.σ[1] = [0.0, 0.0]
m.parameters.γ[1] = [0.0, 0.0]
m.pumptimes = -20.1:19.9

m_data = deepcopy(m)
# scene = plot(m_data)

evaluate!(Z, X, Y, m_data)
z .= reshape(Z, (length(y), length(x)))
z .+= (rand(size(z)...) .- 0.5) ./ 50


m_fitstart = deepcopy(m)
m_fitstart.parameters.A = [4.7, 5.0]
m_fitstart.parameters.ω = [2863.0, 2931.0]
# plot(m_fitstart)

m_fit = deepcopy(m_fitstart)

lower, upper = PumpProbeFitting.generate_bounds(m_fit)

r,g = PumpProbeFitting.fit(z, m_fit, lower, upper)


scene = plot(m_fit, z)

scale!(scene, 1, 1, 400)
update_cam!(scene, FRect3D(Vec3f0(0), Vec3f0(1)))
scene
