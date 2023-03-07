const Ly   = 70
const h    = 1000.0
const ΔB   = 6.0e-2 
const ΔT   = 30.0
const fact = 5.0

"""
    function zonal_wind_stress(y, mid_wind)

returns the zonal wind as per https://egusphere.copernicus.org/preprints/2022/egusphere-2022-186/egusphere-2022-186.pdf
as a function of latitude `y` with  `mid_wind` the wind at the equator (`y = 0.0`)
    
"""
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

"""
    function salinity_flux(y, mid_flux)

returns the salinity flux as a function of latitude `y` 
(similar to https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020gl089135)
"""
@inline function salinity_flux(y)
    if y < -20
        return cubic_profile(y, -70.0, -20.0, -2e-8, 2e-8, 0.0, 0.0) .* 35.0
    elseif y < 0
        return cubic_profile(y, -20.0, 0.0, 2e-8, -4e-8, 0.0, 0.0) .* 35.0
    elseif y < 20
        return cubic_profile(y, 0.0, 20.0, -4e-8, 2e-8, 0.0, 0.0) .* 35.0
    else
        return cubic_profile(y, 20.0, 70.0, 2e-8, -2e-8, 0.0, 0.0) .* 35.0
    end
end

#####
##### Functions specifying initial conditions
#####

@inline exponential_profile(z; Δ = ΔB, Lz = Lz, h = h) = ( Δ * (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) )

@inline parabolic_scaling(y) = - 1 / 70^2 * y^2 + 1
@inline atan_scaling(y)      = (atan(fact*((Ly + y)/Ly - 0.5)) / atan(fact * 0.5) + 1) /2

@inline initial_buoyancy_tangent(x, y, z)  = exponential_profile(z) * atan_scaling(y)
@inline initial_buoyancy_parabola(x, y, z) = exponential_profile(z) * parabolic_scaling(y) 

@inline initial_temperature_parabola(x, y, z) = exponential_profile(z; Δ = ΔT) * parabolic_scaling(y)

@inline function initial_salinity(y, mid_salinity)
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

#####
##### Bottom drag boundary conditions
#####

@inline ϕ²(i, j, k, grid, ϕ) = ϕ[i, j, k]^2

@inline speedᶠᶜᶜ(i, j, k, grid, fields) = (fields.u[i, j, k]^2 + ℑxyᶠᶜᵃ(i, j, k, grid, ϕ², fields.v))^0.5
@inline speedᶜᶠᶜ(i, j, k, grid, fields) = (fields.v[i, j, k]^2 + ℑxyᶜᶠᵃ(i, j, k, grid, ϕ², fields.u))^0.5

@inline u_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, 1] * speedᶠᶜᶜ(i, j, 1, grid, fields)
@inline v_bottom_drag(i, j, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, 1] * speedᶜᶠᶜ(i, j, 1, grid, fields)

@inline u_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.u[i, j, k] * speedᶠᶜᶜ(i, j, k, grid, fields) 
@inline v_immersed_bottom_drag(i, j, k, grid, clock, fields, μ) = @inbounds - μ * fields.v[i, j, k] * speedᶜᶠᶜ(i, j, k, grid, fields) 

#####
##### Wind stress boundary condition
#####

@inline surface_wind_stress(i, j, grid, clock, fields, τ) = τ[j]

#####
##### Tracers boundary condition
#####

@inline function buoyancy_top_relaxation(i, j, grid, clock, fields, p) 

    b = fields.b[i, j, grid.Nz]
    x, y, z = node(Center(), Center(), Center(), i, j, grid.Nz, grid)

    return @inbounds p.λ * (b - p.initial_buoyancy(x, y, z))
end

@inline function temperature_top_relaxation(i, j, grid, clock, fields, p) 

    T  = fields.T[i, j, grid.Nz]
    x, y, z = node(Center(), Center(), Center(), i, j, grid.Nz, grid)
    Trestoring = p.initial_temperature(x, y, z)

    return @inbounds p.λ * (T - Trestoring)
end

@inline function salinity_top_relaxation(i, j, grid, clock, fields, p) 

    S  = fields.S[i, j, grid.Nz]
    Srestoring = p.Ss[j]
    Sflux      = p.Fs[j]

    return @inbounds p.λ * (S - Srestoring) - Sflux
end

#####
##### Utility to concretize a boundary function into an array
#####

@inline function grid_specific_array(boundary_function, grid; scaling = 1000.0)

    Ny   = size(grid, 2)
    arch = architecture(grid)
    
    φ_grid = grid.φᵃᶜᵃ[1:Ny]

    τw = zeros(Ny)
    for (j, φ) in enumerate(φ_grid)
        τw[j] = boundary_function(φ, 0.0) ./ scaling
    end

    return arch_array(arch, -τw)
end