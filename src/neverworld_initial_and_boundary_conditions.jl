const Ly   = 70
const h    = 1000.0
const ΔB   = 6.0e-2 
const ΔT   = 30.0
const fact = 5.0

@inline function zonal_wind_stress(y; τₚ = [0.2, -0.1, -0.02, -0.1, 0.1])
    if y < -45
        return cubic_profile(y, -70.0, -45.0, 0.0, τₚ[1], 0.0, 0.0)
    elseif y < -15
        return cubic_profile(y, -45.0, -15.0, τₚ[1], τₚ[2], 0.0, 0.0)
    elseif y < 0
        return cubic_profile(y, -15.0, 0.0, τₚ[2], τₚ[3], 0.0, 0.0)
    elseif y < 15
        return cubic_profile(y, 0.0, 15.0, τₚ[3], τₚ[4], 0.0, 0.0)
    elseif y < 45
        return cubic_profile(y, 15.0, 45.0, τₚ[4], τₚ[5], 0.0, 0.0)
    else
        return cubic_profile(y, 45.0, 70.0, τₚ[5], 0.0, 0.0, 0.0)
    end
end

@inline exponential_profile(z; Δ = ΔB, Lz = Lz, h = h) = ( Δ * (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) )

@inline parabolic_scaling(y) = - 1 / 70^2 * y^2 + 1
@inline atan_scaling(y)      = (atan(fact*((Ly + y)/Ly - 0.5)) / atan(fact * 0.5) + 1) /2

@inline initial_buoyancy_tangent(x, y, z)  = exponential_profile(z) * atan_scaling(y)
@inline initial_buoyancy_parabola(x, y, z) = exponential_profile(z) * parabolic_scaling(y) 

@inline initial_temperature_parabola(x, y, z) = ΔT * parabolic_scaling(y)

@inline restoring_temperature_parabola(y) = parabolic_scaling(y)

@inline function restoring_salinity_piecewise_cubic(y)
    if y < -20
        return cubic_profile(y, -70.0, -20.0, 34.0, 37.0, 0.0, 0.0)
    elseif y < 0
        return cubic_profile(y, -20.0, 0.0, 37.0, 35.0, 0.0, 0.0)
    elseif y < 20
        return cubic_profile(y, 0.0, 20.0, 35.0, 37.0, 0.0, 0.0)
    else
        return cubic_profile(y, 20.0, 70.0, 37.0, 34.0, 0.0, 0.0)
    end
end

@inline function salinity_flux(y; Sₚ = [-2e-8, 2e-8, -4e-8, 2e-8, -2e-8])
    if y < -20
        return cubic_profile(y, -70.0, -20.0, Sₚ[1], Sₚ[2], 0.0, 0.0) .* 35.0
    elseif y < 0
        return cubic_profile(y, -20.0, 0.0, Sₚ[2], Sₚ[3], 0.0, 0.0) .* 35.0
    elseif y < 20
        return cubic_profile(y, 0.0, 20.0, Sₚ[3], Sₚ[4], 0.0, 0.0) .* 35.0
    else
        return cubic_profile(y, 20.0, 70.0, Sₚ[4], Sₚ[5], 0.0, 0.0) .* 35.0
    end
end

