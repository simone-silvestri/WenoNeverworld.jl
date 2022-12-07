
using Oceananigans.Operators
using Oceananigans.Advection: WENOVectorInvariant, 
		    _left_biased_interpolate_yᵃᶜᵃ,
		    _left_biased_interpolate_xᶜᵃᵃ,
		   _right_biased_interpolate_yᵃᶜᵃ,
		   _right_biased_interpolate_xᶜᵃᵃ,
		    	    upwind_biased_product

function WENO(vector_invariant::SmoothnessStencil, FT::DataType=Float64; 
              order = 5,
              grid = nothing, 
              zweno = true,
              bounds = nothing)

    if !(grid isa Nothing) 
    FT = eltype(grid)
    end

    mod(order, 2) == 0 && throw(ArgumentError("WENO reconstruction scheme is defined only for odd orders"))

    if order < 3
    # WENO(order = 1) is equivalent to UpwindBiased(order = 1)
    return UpwindBiased(order = 1)
    else
    VI = typeof(vector_invariant)
    N  = Int((order + 1) ÷ 2)

    weno_coefficients = compute_reconstruction_coefficients(grid, FT, :WENO; order = N)
    buffer_scheme     = WENO(vector_invariant, FT; grid, order = order - 2, zweno, bounds)
    advecting_velocity_scheme = Centered(FT; grid, order = order - 1)
    end

    return WENO{N, FT, VI, zweno}(weno_coefficients..., bounds, buffer_scheme, advecting_velocity_scheme)
end
