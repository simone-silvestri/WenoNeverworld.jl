using Oceananigans.Operators: ζ₃ᶠᶠᶜ
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Models.HydrostaticFreeSurfaceModels: ZStarSpacingGrid

using Oceananigans.Models: AbstractModel
using Oceananigans.DistributedComputations

const DistributedSimulation = Simulation{<:AbstractModel{<:Distributed}}

maybe_distributed_filename(simulation, output_prefix) = output_prefix
maybe_distributed_filename(sim::DistributedSimulation, output_prefix) = output_prefix * "_$(sim.model.architecture.local_rank)"

"""	
    function standard_outputs!(simulation, output_prefix; overwrite_existing = true, 	
                                                          checkpoint_time    = 100days,	
                                                          snapshot_time      = 30days,	
                                                          surface_time       = 5days,	
                                                          average_time       = 30days,	
                                                          average_window     = average_time,	    
                                                          average_stride     = 10)	

attaches four `JLD2OutputWriter`s to `simulation` with prefix `output_prefix`	

Outputs attached	
================	

- `snapshots` : snapshots of `u`, `v`, `w` and `b` saved every `snapshot_time`	
- `surface_fields` : snapshots of `u`, `v`, `w` and `b` at the surface saved every `surface_time`	
- `averaged_fields` : averages of `u`, `v`, `w`, `b`, `ζ`, `ζ2`, `u2`, `v2`, `w2`, `b2`, `ub`, `vb`, and `wb` 	
                      saved every `average_time` with a window of `average_window` and stride of `average_stride`	
- `checkpointer` : checkpointer saved every `checkpoint_time`	

"""
function standard_outputs!(simulation, output_prefix; overwrite_existing = true, 
                                                      checkpoint_time    = 100days,
                                                      snapshot_time      = 30days,
                                                      surface_time       = 5days,
                                                      average_time       = 30days,
                                                      average_window     = average_time,
                                                      average_stride     = 10)

    output_prefix = maybe_distributed_filename(simulation, output_prefix)
    
    model = simulation.model
    grid  = model.grid

    u, v, w = model.velocities
    b = model.tracers.b

    output_fields = (; u, v, w, b)

    if grid isa ZStarSpacingGrid
        output_fields = merge(output_fields, (; sⁿ = grid.Δzᵃᵃᶠ.sⁿ, wᵍ = grid.Δzᵃᵃᶠ.∂t_∂s))
    end

    u2 = u^2
    v2 = v^2
    b2 = b^2
    w2 = w^2
    vb = v * b
    ub = u * b
    wb = w * b

    ζ  = KernelFunctionOperation{Face, Face, Center}(ζ₃ᶠᶠᶜ, grid, u, v)
    ζ2 = ζ^2

    averaged_fields = (; u, v, w, b, ζ, ζ2, u2, v2, w2, b2, ub, vb, wb)

    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                                  schedule = TimeInterval(snapshot_time),
                                                                  filename = output_prefix * "_snapshots",
                                                                  overwrite_existing)

    simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, output_fields;
                                                                  schedule = TimeInterval(surface_time),
                                                                  filename = output_prefix * "_surface",
                                                                  indices = (:, :, grid.Nz),
                                                                  overwrite_existing)

    simulation.output_writers[:averaged_fields] = JLD2OutputWriter(model, averaged_fields;
                                                                   schedule = AveragedTimeInterval(average_time, window=average_window, stride = average_stride),
                                                                   filename = output_prefix * "_averages",
                                                                   overwrite_existing)

    simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                            schedule = TimeInterval(checkpoint_time),
                                                            prefix = output_prefix * "_checkpoint",
                                                            overwrite_existing)

    return nothing
end


"""	
    function checkpoint_outputs!(simulation, output_prefix; overwrite_existing = true, checkpoint_time = 100days)	
        
attaches a `Checkpointer` to the simulation with prefix `output_prefix` that is saved every `checkpoint_time`	
"""
function checkpoint_outputs!(simulation, output_prefix; overwrite_existing = true, checkpoint_time = 100days)

    output_prefix = maybe_distributed_filename(simulation, output_prefix)

    model = simulation.model

    simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                            schedule = TimeInterval(checkpoint_time),
                                                            prefix = output_prefix * "_checkpoint",
                                                            overwrite_existing)

    return nothing
end

"""	
    function reduced_outputs!(simulation, output_prefix; overwrite_existing = true, 	
                                                         checkpoint_time    = 100days,	
                                                         snapshot_time      = 30days,	
                                                         surface_time       = 1days,	
                                                         bottom_time        = 1days)	

attaches four `JLD2OutputWriter`s to `simulation` with prefix `output_prefix`	

Outputs attached	
================	

- `snapshots` : snapshots of `u`, `v`, `w` and `b` saved every `snapshot_time`	
- `surface_fields` : snapshots of `u`, `v`, `w` and `b` at the surface saved every `surface_time`	
- `bottom_fields` : snapshots of `u`, `v`, `w` and `b` at the bottom (`k = 2`) saved every `bottom_time`	
- `checkpointer` : checkpointer saved every `checkpoint_time`	

"""
function reduced_outputs!(simulation, output_prefix; overwrite_existing = true, 
                                                     checkpoint_time    = 100days,
                                                     snapshot_time      = 30days,
                                                     surface_time       = 1days,
                                                     bottom_time        = 1days)

    output_prefix = maybe_distributed_filename(simulation, output_prefix)

    model = simulation.model
    grid  = model.grid

    u, v, w = model.velocities
    b = model.tracers.b

    output_fields = (; u, v, w, b)

    if grid isa ZStarSpacingGrid
        output_fields = merge(output_fields, (; sⁿ = grid.Δzᵃᵃᶠ.sⁿ, wᵍ = grid.Δzᵃᵃᶠ.∂t_∂s))
    end
    
    simulation.output_writers[:snapshots] = JLD2OutputWriter(model, output_fields;
                                                                schedule = TimeInterval(snapshot_time),
                                                                filename = output_prefix * "_snapshots",
                                                                overwrite_existing)
   

    simulation.output_writers[:surface_fields] = JLD2OutputWriter(model, output_fields;
                                                                    schedule = TimeInterval(surface_time),
                                                                    filename = output_prefix * "_surface",
                                                                    indices = (:, :, grid.Nz),
                                                                    overwrite_existing)
                                                                                                                                
    simulation.output_writers[:bottom_fields] = JLD2OutputWriter(model, output_fields;
                                                                    schedule = TimeInterval(bottom_time),
                                                                    filename = output_prefix * "_bottom",
                                                                    indices = (:, :, 2),
                                                                    overwrite_existing)
    
    simulation.output_writers[:checkpointer] = Checkpointer(model;
                                                            schedule = TimeInterval(checkpoint_time),
                                                            prefix = output_prefix * "_checkpoint",
                                                            overwrite_existing)

end                                                 

"""	
    function vertically_averaged_outputs!(simulation, output_prefix; overwrite_existing = true, 	
                                                                     average_time       = 30days,	
                                                                     average_window     = average_time,	    
                                                                     average_stride     = 10)	

attaches a `JLD2OutputWriter`s to `simulation` with prefix `output_prefix`	

Outputs attached	
================	

- `vertically_averaged_outputs` : average of `KE` and heat content (integral of temperature in ᵒCm³)

"""
function vertically_averaged_outputs!(simulation, output_prefix; overwrite_existing = false, 
                                                                 average_time       = 30days,
                                                                 average_window     = average_time,
                                                                 average_stride     = 10)

    output_prefix = maybe_distributed_filename(simulation, output_prefix)

    model = simulation.model
    g = simulation.model.free_surface.gravitational_acceleration
    α = 2e-4
    u, v, _ = model.velocities
    
    T  = Field(g * α * model.tracers.b)
    KE = Field(0.5 * (u^2 + v^2))

    tke_average  = Average(KE,  dims = 3)
    heat_content = Integral(T,  dims = 3)

    output_fields = (; tke_average, heat_content)
                                                      
    simulation.output_writers[:vertically_averaged_outputs] = JLD2OutputWriter(model, output_fields;
                                                              schedule = AveragedTimeInterval(average_time, window=average_window, stride = average_stride),
                                                              filename = output_prefix * "_vertical_average",
                                                              overwrite_existing)
end