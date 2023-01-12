const Ly   = 70
const h    = 1000.0
const ΔB   = 6.0e-2 
const fact = 5.0

@inline function zonal_wind_stress(y, mid_wind)
    if y < -45
        return cubic_profile(y, -70.0, -45.0, 0.0, 0.2, 0.0, 0.0)
    elseif y < -15
        return cubic_profile(y, -45.0, -15.0, 0.2, -0.1, 0.0, 0.0)
    elseif y < 0
        return cubic_profile(y, -15.0, 0.0, -0.1, mid_wind, 0.0, 0.0)
    elseif y < 15
        return cubic_profile(y, 0.0, 15.0, mid_wind, -0.1, 0.0, 0.0)
    elseif y < 45
        return cubic_profile(y, 15.0, 45.0, -0.1, 0.1, 0.0, 0.0)
    else
        return cubic_profile(y, 45.0, 70.0, 0.1, 0.0, 0.0, 0.0)
    end
end

@inline exponential_profile(z; Δ = ΔB, Lz = Lz, h = h) = ( Δ * (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) )

@inline parabolic_scaling(y) = - 1 / 70^2 * y^2 + 1
@inline atan_scaling(y)      = (atan(fact*((Ly + y)/Ly - 0.5)) / atan(fact * 0.5) + 1) /2

@inline initial_buoyancy_tangent(x, y, z)  = exponential_profile(z) * atan_scaling(y)
@inline initial_buoyancy_parabola(x, y, z) = exponential_profile(z) * parabolic_scaling(y) 

const σᵤ = 0.5
const σᵥ = 50.0

const x₁ =   30
const y₁ = - 50
const z₁ = - 500

@inline initial_tracer(x, y, z) = exp( - ((x - x₁)^2 + (y - y₁)^2) / (2*σᵤ^2) - (z - z₁)^2 / (2 * σᵥ^2))

const i₀ = 130
const j₀ = 85
const k₀ = 25

const σᵢₕ = 10
const σᵢᵥ = 5

@inline initial_tracer_idx(i, j, k) = exp( - ((i - i₀)^2 + (j - j₀)^2) / (2*σᵢₕ^2) - (k - k₀)^2 / (2 * σᵢᵥ^2))

using Oceananigans.Utils: launch!
using Oceananigans.Grids: architecture

@kernel function _set_c_field!(c)
    i, j, k = @index(Global, NTuple)
    @inbounds c[i, j, k] = initial_tracer_idx(i, j, k)
end

@inline function set_c_field!(c)
    arch  = architecture(c)
    grid  = c.grid
    event = launch!(arch, grid, :xyz, _set_c_field!, c)

    wait(event)
end
