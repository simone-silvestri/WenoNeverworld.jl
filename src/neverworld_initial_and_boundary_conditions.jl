const Lz   = 4000.0
const Ly   = 70
const h    = 1000.0
const ΔB   = 6.0e-2 
const fact = 5.0

@inline function wind_stress(y, mid_wind)
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

@inline atan_scaling(y)           = (atan(fact*((Ly + y)/Ly - 0.5)) / atan(fact * 0.5) + 1) /2
@inline initial_buoyancy(x, y, z) = ( ΔB * (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) ) * atan_scaling(y)