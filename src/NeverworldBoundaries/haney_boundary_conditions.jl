using WenoNeverworld.Auxiliaries: parabolic_scaling
using Oceananigans.Grids: node
using Base

####
#### Some utilities to apply the BC to the correct tracer
####

struct Buoyancy end
struct Temperature end
struct Salinity end

@inline zerofunc(args...) = 0

# Only tracers for the moment: all at `Center`s
@inline field_location(varname) = (Center(), Center(), Center())

@propagate_inbounds Base.getindex(fields, i, j, k, ::Buoyancy)    = @inbounds fields.b[i, j, k]
@propagate_inbounds Base.getindex(fields, i, j, k, ::Temperature) = @inbounds fields.T[i, j, k]
@propagate_inbounds Base.getindex(fields, i, j, k, ::Salinity)    = @inbounds fields.S[i, j, k]

@inline getvalue(a::Number,        args...) = a
@inline getvalue(a::AbstractArray, i, j, k, args...) = @inbounds a[i, j, k]

@inline function getvalue(a::Function, i, j, k, grid, loc, time, args...) 
    X = node(i, j, k, grid, loc...)
    return a(X..., time, args...)
end

struct HaneyBoundaryCondition{F, T, B, V, P} <: Function
    flux :: F
    pumping_velocity :: T
    restoring_profile :: B
    varname :: V
    parameters :: P
end

"""
    HaneyBoundaryCondition(; flux = zerofunc,
                             pumping_velocity = 0.0,
                             restoring_profile = zerofunc,
                             varname = Temperature(),
                             parameters = nothing)

Constructs a HaneyBoundaryCondition object, which includes a prescribed `flux`
and a restoring to a `restoring_profile` with a `pumping_velocity`

## Arguments
- `flux`: The flux function.
- `pumping_velocity`: The pumping velocity.
- `restoring_profile`: The restoring profile function.
- `varname`: The variable name.
- `parameters`: Additional parameters.

## Returns
A HaneyBoundaryCondition object.
"""
function HaneyBoundaryCondition(; flux = zerofunc,
                                  pumping_velocity = 0.0,
                                  restoring_profile = zerofunc,
                                  varname = Temperature(),
                                  parameters = nothing)
                                  
    return HaneyBoundaryCondition(flux, pumping_velocity, restoring_profile, varname, parameters)
end

# The function called by `apply_bc!`
function (bc::HaneyBoundaryCondition)(i, j, grid, clock, fields)

    kᴺ = grid.Nz

    # First retrieve the prescribed flux
    loc  = field_location(bc.varname)
    flux = getvalue(bc.flux, i, j, kᴺ, grid, loc, clock.time, bc.parameters)

    # Now we calculate the restoring
    var  = @inbounds fields[i, j, kᴺ, bc.varname] 
    var★ = getvalue(bc.restoring_profile, i, j, kᴺ, grid, loc, clock.time, bc.parameters)

    return flux + bc.pumping_velocity * (var - var★)
end

Adapt.adapt_structure(to, b::HaneyBoundaryCondition) = 
    HaneyBoundaryCondition(Adapt.adapt(to, b.flux),
                           Adapt.adapt(to, b.pumping_velocity),
                           Adapt.adapt(to, b.restoring_profile),
                           Adapt.adapt(to, b.varname),
                           Adapt.adapt(to, b.params))
