module PumpProbeFitting

using PumpProbeModels
using LsqFit
using Printf

mutable struct PumpProbeFit
    m::PumpProbeModel
    f # Fit
    d # Data
end

function fit(z, m::PumpProbeModel, lower::PumpProbeModels.ParamContainer,
                                   upper::PumpProbeModels.ParamContainer)

    # Wrapper for the cuve_fit function
    function model(F, xt, p)
        pack_parameters!(m, p)
        x = @view xt[:,1]
        t = @view xt[:,2]
        evaluate!(F, x, t, m)
        F
    end

    # convert wavenumber and delaytime inputs to 1D
    x = eltype(m.wavenumbers)[]
    t = eltype(m.delaytimes)[]
    for _x in m.wavenumbers, _t in m.delaytimes
        push!(x, _x)
        push!(t, _t)
    end
    xt = [x t]
    Z = z[:]

    start = deepcopy(m.parameters)
    p0 = unpack_parameters(start)
    pl = unpack_parameters(lower)
    pu = unpack_parameters(upper)

    f = curve_fit(model, xt, Z, p0, lower=pl, upper=pu, maxIter=500,
                  autodiff=:finiteforward, inplace=true, show_trace=false)

    pack_parameters!(m, f.param)
    fitparams = m.parameters

    printparams(start, lower, upper, fitparams)

    PumpProbeFit(m, f, z), model(z[:], xt, f.param)
end

function printparams(s, l, u, f)
    param_names = fieldnames(typeof(s))
    for n in param_names
        @printf "%s\n" n
        if isa(getfield(s, n), AbstractArray)
            for i = 1:length(getfield(s, n))
                if isa(getfield(s, n)[i], AbstractArray)
                    for j = 1:length(getfield(s, n)[i])
                        @printf "%8.3f\t%8.3f\t%8.3f\t%8.3f\n" getfield(s, n)[i][j] getfield(l, n)[i][j] getfield(u, n)[i][j] getfield(f, n)[i][j]
                    end
                else
                    @printf "%8.3f\t%8.3f\t%8.3f\t%8.3f\n" getfield(s, n)[i] getfield(l, n)[i] getfield(u, n)[i] getfield(f, n)[i]
                end
            end
        else
            @printf "%8.3f\t%8.3f\t%8.3f\t%8.3f\n" getfield(s, n) getfield(l, n) getfield(u, n) getfield(f, n)
        end
        @printf "________________________________________________\n"
    end
end

function generate_bounds(m::PumpProbeModel)
    l = deepcopy(m.parameters)
    u = deepcopy(m.parameters)

    l.A    .= 1.0
    l.a    .= -0.5
    l.τ[1] .= 30.0
    l.ω   .*= 0.9
    l.Δω   .= 0.0
    l.σ[1] .= 0.0
    l.Γ    .= 5.0
    l.ΔΓ   .= 0.0
    l.γ[1] .= 0.0
    l.t0   -= 30.0

    u.A    .= 20.0
    u.a    .= 1.0
    u.τ[1] .= 500.0
    u.ω   .*= 1.1
    u.Δω   .= 0.0
    u.σ[1] .= 0.0
    u.Γ    .= 15.0
    u.ΔΓ   .= 0.0
    u.γ[1] .= 0.0
    u.t0   += 30.0

    return l, u
end


unpack_parameters(m::PumpProbeModel) = unpack_parameters(m.parameters)

function unpack_parameters(p::PumpProbeModels.ParamContainer)
    if length(p.τ) == 1 # :exp decay mode
        params_flat = [
            p.A...,
            p.a...,
            p.τ[1]...,
            p.ω...,
            p.Δω...,
            p.σ[1]...,
            p.Γ...,
            p.ΔΓ...,
            p.γ[1]...,
            p.t0
        ]
    elseif length(p.τ) == 2 # :doubleexp decay mode
        params_flat = [
            p.A...,
            p.a...,
            p.τ[1]...,
            p.τ[2]...,
            p.ω...,
            p.Δω...,
            p.σ[1]...,
            p.Γ...,
            p.ΔΓ...,
            p.γ[1]...,
            p.t0
        ]
    end
    return params_flat
end


function pack_parameters!(m::PumpProbeModel, fp)
    p = m.parameters
    if m.decaymode == :exp
        # Get number of resonances
        n = (length(fp) - 1) ÷ 9
        p.A    .= fp[   1: n]
        p.a    .= fp[ n+1:2n]
        p.τ[1] .= fp[2n+1:3n]
        p.ω    .= fp[3n+1:4n]
        p.Δω   .= fp[4n+1:5n]
        p.σ[1] .= fp[5n+1:6n]
        p.Γ    .= fp[6n+1:7n]
        p.ΔΓ   .= fp[7n+1:8n]
        p.γ[1] .= fp[8n+1:9n]
        p.t0    = fp[9n+1]
    elseif m.decaymode == :doubleexp
        # TODO: Implement fully
        error("not implemented")
        params_flat = [
            p.A...,
            p.a...,
            p.τ[1]...,
            p.τ[2]...,
            p.ω...,
            p.Δω...,
            p.σ[1]...,
            p.Γ...,
            p.ΔΓ...,
            p.γ[1]...,
            p.t0
        ]
    end
    return m
end

end # module
