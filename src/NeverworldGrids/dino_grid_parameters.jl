#####
##### Functions that build a grid equivalent to the one used in the DINO simulation
#####

# Code reference: https://github.com/vopikamm/DINO

# Convenience function to pass DINO parameters to the 
# grid constructor as `NeverworldGrid(res; dino_parameters(res)...)`
function dino_parameters(resolution; longitude = (0, 50))
    latitude   = dino_latitude(resolution)
    z_faces    = dino_vertical_coordinate()
    bathymetry = DinoBathymetry(; resolution)

    return (; longitude, latitude, z_faces, bathymetry)
end

"""
    dino_latitude(resolution; φ_max = 70)

a stretched latitudinal coordinate, following the mercator projections,
with smaller Δφ at the poles
"""
function dino_latitude(resolution; φ_max = 70) 
    Δλ = resolution

    φ₀ = 0
    φ₊ = eltype(Δλ)[]
    j  = 1
    φ₀ = 180 / π * asin(tanh(π / 180 * Δλ * j))
    while φ₀ <= φ_max
        push!(φ₊, φ₀)
        j += 1
        φ₀ = 180 / π * asin(tanh(π / 180 * Δλ * j))
    end

    φ₋ = - reverse(φ₊)
    φ  = vcat(φ₋, zero(Δλ), φ₊)

    return φ
end
       
# A stretched vertical coordinate. See https://github.com/vopikamm/DINO
function dino_vertical_coordinate(; Nz     = 36,     # number of vertical layers
                                    Lz     = 4000.0, # depth of the domain
                                    kᵗʰ    = 35,     
                                    aᶜʳ    = 10.5,
                                    Δz_min = 10)

    K = Nz+1
                            
    z_faces = zeros(Nz+1)
    
    a₀ = (Δz_min - Lz / (K - 1)) / (
          tanh((1 - kᵗʰ) / aᶜʳ) -
          aᶜʳ / (K - 1) * (log(cosh((K - kᵗʰ) / aᶜʳ)) -
                           log(cosh((1  - kᵗʰ) / aᶜʳ)))
    )

    a₁ = Δz_min - a₀ * tanh((1 - kᵗʰ) / aᶜʳ)
    a₂ = -a₁ - a₀ * aᶜʳ * log(cosh((1 - kᵗʰ) / aᶜʳ))
    
    for k in 1:Nz+1
        z_faces[k] = a₂ + a₁ * k + a₀ * aᶜʳ * log(cosh((k - kᵗʰ) / aᶜʳ)) 
    end

    return reverse( - z_faces)
end

Base.@kwdef struct DinoBathymetry
    λ_minimum :: Float64       = 0
    λ_maximum :: Float64       = 50
    φ_minimum :: Float64       = -70
    φ_maximum :: Float64       = 70
    resolution :: Float64      = 1/4
    φ_channel_min :: Float64   = -65   # Minimum channel latitude on tracer-point (approx.)     
    φ_channel_max :: Float64   = -45   # Maximum channel latitude on tracer-point (approx.)
    slope :: Float64           = 3     # slope around the vertical walls
    slope_sill :: Float64      = 4     # slope around the Scotia arc
    H_max :: Float64           = 4000  # Maximum depth of the bathymetry on w-velocity-point
    H_min :: Float64           = 2000  # Minimum depth of the bathymetry on w-velocity-point (approx.)
    H_sill ::Float64           = 2500  # Depth of the Scotia arc sill [meters]
end

function (params::DinoBathymetry)(λ, φ, args...)
    λ⁻  = params.λ_minimum
    λ⁺  = params.λ_maximum
    φ⁻  = params.φ_minimum
    φ⁺  = params.φ_maximum
    φᶜ⁻ = params.φ_channel_min
    φᶜ⁺ = params.φ_channel_max 
    H⁺  = params.H_max
    H⁻  = params.H_min
    Hs  = params.H_sill
    𝒮   = params.slope
    𝒮s  = params.slope_sill
    Δλ  = params.resolution

    Δλᵂ = abs(λ⁺ - λ⁻)
    Δφᶜ = abs(φᶜ⁺ - φᶜ⁻)

    # Start by applying the slopes around the channel 
    # and around the solid walls
    channel   = (1 - exp_bathymetry(φ, φᶜ⁻, φᶜ⁺, Δλᵂ, 𝒮, Δφᶜ / 2))
    xboundary = (1 - exp_bathymetry(λ, λ⁻,  λ⁺,  Δλᵂ, 𝒮, Δφᶜ / 2))
    outside_channel = (φ ≤ φᶜ⁻) | (φ ≥ φᶜ⁺)
    zx = 1 - (xboundary * channel + outside_channel * xboundary)

    𝒮l = cos(π * φ⁺ /180) * 𝒮
    
    zy = exp_bathymetry(φ, φ⁻, φ⁺, Δλᵂ, 𝒮l, Δφᶜ / 2)
    bathymetry = zx * zy * (H⁺ - H⁻) + H⁻

    # Add the gaussian ring that represents the Scotia Arch
    φ₀ = (φᶜ⁺ + φᶜ⁻) / 2
    radius = abs(φᶜ⁺ - φᶜ⁻) / 2

    # taper the gaussian ring eastward 
    taper = smooth_step(λ, λ⁻, λ⁻ + 𝒮s)
    
    # shift grid
    zx = λ - λ⁻
    zy = φ - φ₀

    # Calculate the exponent part of the Gaussian function
    exp_arg = (- zx^2 - zy^2 + 2 * radius * sqrt(zx^2 + zy^2) - radius^2) / 𝒮s^2

    # Calculate the resulting bathymetry depth
    gauss_ring = ifelse(bathymetry ≥ Hs, 
                        (Hs - bathymetry) * exp(exp_arg) + bathymetry, 
                        bathymetry)
    
    # update bathymetry
    bathymetry = taper * gauss_ring + (1 - taper) * bathymetry

    # Leave space for a channel but fill in all at least one longitude point with a vertical wall
    inside_wall = ((λ ≤ λ⁻ + Δλ / 2) | (λ ≥ λ⁺ - Δλ / 2)) & ((φ < φᶜ⁻) | (φ > φᶜ⁺))
    bathymetry  = ifelse(inside_wall, zero(λ), bathymetry)

    # Convert to a depth (negative values)
    return - bathymetry
end

# A smoothing function that transitions between
# 0 where m < mᴸ to 1 where m > mᴿ
function smooth_step(m, mᴸ, mᴿ)
    x    = (m - mᴸ) / (mᴿ - mᴸ)
    step = 6 * x^5 - 15 * x^4 + 10 * x^3 
    return ifelse(m < mᴸ, zero(m),
           ifelse(m > mᴿ, one(m), step))
end

# Exponential function to calculate smooth slopes around vertical walls.
function exp_bathymetry(m, mᴸ, mᴿ, Δλ, 𝒮, Δm)

    n = 1 + exp(- Δλ / 𝒮)

    taper_left  = (m ≥ mᴸ)      & (m ≤ mᴸ + Δm)
    taper_right = (m ≥ mᴿ - Δm) & (m ≤ mᴿ)

    step_left  =     smooth_step(m, mᴸ, mᴸ + Δm)
    step_right = 1 - smooth_step(m, mᴿ - Δm, mᴿ)

    bat_left  = (1 - exp(-(m - mᴸ) / 𝒮) / n) * (1 - step_left)  + step_left
    bat_right = (1 - exp( (m - mᴿ) / 𝒮) / n) * (1 - step_right) + step_right

    bat = ifelse(taper_right, bat_right, 
          ifelse(taper_left,  bat_left, one(n)))
    
    return bat
end
