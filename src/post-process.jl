import Pkg
Pkg.activate("/home/ssilvest/stable_oceananigans/Oceananigans.jl/")

using CairoMakie
using JLD2
using Oceananigans
using Oceananigans.Grids: on_architecture
using Oceananigans.AbstractOperations: Az, GridMetricOperation
using Oceananigans.Operators: ζ₃ᶠᶠᶜ
using Oceananigans.Coriolis: fᶠᶠᵃ
using Oceananigans.Models.HydrostaticFreeSurfaceModels: VerticalVorticityField
include("initial_conditions.jl")

grid = on_architecture(CPU(), grid)

file = jldopen("test.jld2")
iterations = parse.(Int, keys(file["timeseries/t"]))

coriolis = HydrostaticSphericalCoriolis()

u  = Vector(undef, length(iterations))
v  = Vector(undef, length(iterations))
δ  = Vector(undef, length(iterations))
ζ  = Vector(undef, length(iterations))
η  = Vector(undef, length(iterations))
ν  = Vector(undef, length(iterations))
PV = Vector(undef, length(iterations))

Cν = 0.2 / (8π^2)

@inline PV_func(i, j, k, grid, coriolis, u, v) = fᶠᶠᵃ(i, j, k, grid.underlying_grid, coriolis) + ζ₃ᶠᶠᶜ(i, j, k, grid.underlying_grid, u, v)

for i in 1:length(iterations)
    @info "Reading iteration $i of $(length(iterations))"
    u_tmp = file["timeseries/u/" * string(iterations[i])]
    v_tmp = file["timeseries/v/" * string(iterations[i])]
    η_tmp = file["timeseries/η/" * string(iterations[i])]

    u[i] = Field((Face, Center, Nothing), grid)
    v[i] = Field((Center, Face, Nothing), grid)
    η[i] = Field((Center, Center, Nothing), grid)
    
    set!(u[i], u_tmp)
    set!(v[i], v_tmp)
    set!(η[i], η_tmp)

    δ[i] = Field(∂x(u[i]) + ∂y(v[i]))
    PV_op = KernelFunctionOperation{Face, Face, Nothing}(PV_func, grid, computed_dependencies=(coriolis, u[i], v[i]))
    PV[i] = Field(PV_op)

    ζ_op = KernelFunctionOperation{Face, Face, Nothing}(ζ₃ᶠᶠᶜ, grid, computed_dependencies=(u[i], v[i]))
    ζ[i] = Field(ζ_op)

    area = GridMetricOperation((Center, Center, Center), Az, grid)

    ν[i] = Field(Cν * area^2 * (δ[i]^2 + ζ[i]^2)^(0.5))

    compute!(δ[i])
    compute!(ζ[i])
    compute!(ν[i])
end

using Oceananigans.Operators: Δx, Δy
using Oceananigans.Units
@inline νm(i, j, k, grid) = (1 / (1 / Δx(i, j, k, grid, Center(), Center(), Center())^2 + 1 / Δy(i, j, k, grid, Center(), Center(), Center())^2 ))^2 / (1days)

ν_op = KernelFunctionOperation{Center, Center, Nothing}(νm, grid)
νb   = Field(ν_op)
compute!(νb)

iter = Observable(1)

δ_observable = @lift interior(δ[$iter], :, :, 1)
u_observable = @lift interior(u[$iter], :, :, 1)
ζ_observable = @lift interior(ζ[$iter], :, :, 1)
v_observable = @lift interior(v[$iter], :, :, 1)
η_observable = @lift interior(η[$iter], :, :, 1)

fig = Figure(resolution = (1000, 1000))

ax = Axis(fig[1, 1])
hm = heatmap!(ax, u_observable, colorrange = (-0.1, 0.1), colormap = :balance) 
ax = Axis(fig[1, 2])
hm = heatmap!(ax, ζ_observable, colorrange = (-1.0e-5, 1.0e-5), colormap = :blues) 
ax = Axis(fig[2, 1])
hm = heatmap!(ax, v_observable, colorrange = (-0.1, 0.1), colormap = :balance) 
ax = Axis(fig[2, 2])
hm = heatmap!(ax, η_observable, colorrange = (-1, 1), colormap = :viridis) 

record(fig, "surface.mp4", collect(1:length(u)), framerate=12) do i
    @info "Plotting iteration $i of $(length(u))..."
    iter[] = i
end