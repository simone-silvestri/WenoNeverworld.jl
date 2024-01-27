@inline zerofunc(args...) = 0

function validate_tracer_boundary_conditions(tracers, tracer_boundary_conditions)
    for tracer in tracers
        if !(hasproperty(tracer_boundary_conditions, tracer))
            tracer_boundary_conditions = merge(tracer_boundary_conditions, (; tracer => zerofunc))
        end
    end
    return tracer_boundary_conditions
end

materialize_tracer_boundary_conditions(tracers::NamedTuple{(), Tuple{}}, args...) = NamedTuple() 

function materialize_tracer_boundary_conditions(tracers, grid, tracer_bcs)
    bcs = NamedTuple()
    for t in tracers
        bc = getproperty(tracer_bcs, t)
        bc = regularize_boundary_condition(bc, grid)
        top_bc = FluxBoundaryCondition(bc, discrete_form=true)
        bcs = merge(bcs, (; t => FieldBoundaryConditions(top = top_bc)))
    end

    return bcs
end