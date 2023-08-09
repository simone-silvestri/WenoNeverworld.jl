using Oceananigans
using Oceananigans.Units
using WenoNeverworld
using Oceananigans: prognostic_fields
using Oceananigans.BuoyancyModels: g_Earth
using Oceananigans.Grids: minimum_xspacing, minimum_yspacing
using Oceananigans.TurbulenceClosures
using Oceananigans.TurbulenceClosures: ExplicitTimeDiscretization
using Oceananigans.Coriolis: ActiveCellEnstrophyConservingScheme
using JLD2

import Oceananigans.OutputWriters: set!, set_time_stepper!, set_time_stepper_tendencies!

function set!(model, filepath::AbstractString)

    jldopen(filepath, "r") do file

        # Do NOT validate the grid!
        model_fields = prognostic_fields(model)

        for name in (:u, :v, :b, :η)
            if string(name) ∈ keys(file) # Test if variable exist in checkpoint.
                model_field = model_fields[name]
                parent_data = file["$name/data"] #  Allow different halo size by loading only the interior
        		copyto!(model_field.data.parent, parent_data)
            end
        end

        set_time_stepper!(model.timestepper, file, model_fields)

        checkpointed_clock = file["clock"]

        # Update model clock
        model.clock.iteration = checkpointed_clock.iteration
        model.clock.time = checkpointed_clock.time
    end

    return nothing
end

function set_time_stepper_tendencies!(timestepper, file, model_fields)
    for name in (:u, :v, :b)
        if string(name) ∈ keys(file["timestepper/Gⁿ"]) # Test if variable tendencies exist in checkpoint
            # Tendency "n"
            parent_data = file["timestepper/Gⁿ/$name/data"]

            tendencyⁿ_field = timestepper.Gⁿ[name]
            copyto!(tendencyⁿ_field.data.parent, parent_data)

            # Tendency "n-1"
            parent_data = file["timestepper/G⁻/$name/data"]

            tendency⁻_field = timestepper.G⁻[name]
            copyto!(tendency⁻_field.data.parent, parent_data)
        else
            @warn "Tendencies for $name do not exist in checkpoint and could not be restored."
        end
    end

    return nothing
end

output_dir    = joinpath(@__DIR__, "./")
output_dir = "/storage2/"
@show output_prefix = output_dir * "WenoNeverworldData/weno_one"

arch = GPU()
new_degree = 1
old_degree = 1

grid = NeverworldGrid(arch, new_degree, latitude = (-70, 70), H=7)
orig_grid = NeverworldGrid(arch, old_degree, latitude = (-70, 70), H=7)

# Extend the vertical advection scheme
interp_init = true
#init_file = false # "/orcd/nese/raffaele/001/sandre/WenoNeverworld/weno_eighth_checkpoint_iteration7938374.jld2"
init_file = "/storage2/WenoNeverworldData/weno_one_checkpoint_iteration26287125.jld2"
# Simulation parameters
Δt        = 20minutes
stop_time = 3000years

tracer_advection      = WENO(grid.underlying_grid)
momentum_advection    = VectorInvariant(vorticity_scheme = WENO(order = 9), 
                                         vertical_scheme = WENO(grid.underlying_grid))

biharmonic_viscosity  = nothing
convective_adjustment  = RiBasedVerticalDiffusivity()

free_surface = SplitExplicitFreeSurface(; cfl = 0.7, grid)
coriolis = HydrostaticSphericalCoriolis(scheme = ActiveCellEnstrophyConservingScheme())
# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file, 
                                          tracer_advection, momentum_advection, 
                                          biharmonic_viscosity, free_surface, orig_grid, coriolis)

# Adaptable time step
wizard = TimeStepWizard(; cfl = 0.35, max_Δt = 40minutes, min_Δt = 15minutes, max_change = 1.1)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Let's goo!
@info "Running with Δt = $(prettytime(simulation.Δt))"

overwrite_existing = true
# standard_outputs!(simulation, output_prefix; overwrite_existing)
checkpoint_outputs!(simulation, output_prefix; overwrite_existing = false, checkpoint_time = years)

# initializing the time for wall_time calculation
run_simulation!(simulation; interp_init, init_file)

