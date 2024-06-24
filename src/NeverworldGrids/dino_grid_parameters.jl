
function dino_parameters(resolution)
    latitude = dino_latitude(resolution)
    bathymetry = DinoBathymetry()
    longitude = dino_longitude
    z_faces = dino_vertical_coordinate()

    return (; longitude, latitude, z_faces, bathymetry)
end

dino_longitude = (0, 50)

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

    φ = vcat(φ₋, eltype(Δλ)(0), φ₊)

    return φ
end
       
function dino_vertical_coordinate(; Nz     = 36, 
                                    Lz     = 4000.0,
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
    φ_channel_min :: Float64   = -65   # Minimum channel latitude on tracer-point (approx.)     
    φ_channel_max :: Float64   = -45   # Maximum channel latitude on tracer-point (approx.)
    slope :: Float64           = 3     # slope
    H_max :: Float64           = 4000  # Maximum depth of the bathymetry on w-velocity-point
    H_min :: Float64           = 2000  # Minimum depth of the bathymetry on w-velocity-point (approx.)
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
    𝒮   = params.slope

    Δλᵂ = abs(λ⁺ - λ⁻)
    Δφᶜ = abs(φᶜ⁺ - φᶜ⁻)
    
    zy_cha = exp_bathymetry(φ, φᶜ⁻, φᶜ⁺, Δλᵂ, 𝒮, Δφᶜ / 2)
    zx_raw = exp_bathymetry(λ, λ⁻,  λ⁺,  Δλᵂ, 𝒮, Δφᶜ / 2)

    zx = zx_raw * (1 - zy_cha) + zy_cha

    𝒮l = cos(π * φ⁺ /180) * 𝒮
    
    zy = exp_bathymetry(φ, φ⁻, φ⁺, Δλᵂ, 𝒮l, Δφᶜ / 2)
    depth =  - zx * zy * (H⁺ - H⁻) - H⁻

    return depth
end

function exp_bathymetry(m, mᴸ, mᴿ, Δλ, 𝒮, Δm)

    n = 1 + exp(- Δλ / 𝒮)

    taper_left  = (m ≥ mᴸ) & (m ≤ mᴸ + Δm)
    taper_right = (m ≥ mᴿ) & (m ≤ mᴿ + Δm)

    step_left  =     smooth_step(m, mᴸ, mᴸ + Δm)
    step_right = 1 - smooth_step(m, mᴿ - Δm, mᴿ)

    bat_left  = (1 - exp(-(m  - mᴸ) / 𝒮) / n) * (1 - step_left)  + step_left
    bat_right = (1 - exp(-(mᴿ - m)  / 𝒮) / n) * (1 - step_right) + step_right

    bat = ifelse(taper_right, bat_right, 
          ifelse(taper_left,  bat_left, one(n)))

    return bat
end

function smooth_step(m, mᴸ, mᴿ)
    x    = (m - mᴸ) / (mᴿ - mᴸ)
    step = 6 * x^5 - 15 * x^4 + 10 * x^3 
    return ifelse(m < mᴸ, zero(m),
           ifelse(m > mᴿ, one(m), step))
end