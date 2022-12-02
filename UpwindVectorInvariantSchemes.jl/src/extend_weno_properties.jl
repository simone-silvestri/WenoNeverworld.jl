
using Oceananigans.Operators
using Oceananigans.Advection: WENOVectorInvariant, 
		    _left_biased_interpolate_yᵃᶜᵃ,
		    _left_biased_interpolate_xᶜᵃᵃ,
		   _right_biased_interpolate_yᵃᶜᵃ,
		   _right_biased_interpolate_xᶜᵃᵃ,
		    	    upwind_biased_product

using Oceananigans.Advection: SmoothnessStencil


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

# Interpolation functions
for (interp, dir, val, cT) in zip([:xᶠᵃᵃ, :yᵃᶠᵃ, :zᵃᵃᶠ], [:x, :y, :z], [1, 2, 3], [:XT, :YT, :ZT]) 
    for side in (:left, :right)
        interpolate_func = Symbol(:inner_, side, :_biased_interpolate_, interp)
        stencil          = Symbol(side, :_stencil_, dir)
        weno_weights     = Symbol(side, :_biased_weno_weights)
        biased_p         = Symbol(side, :_biased_p)
        
        @eval begin
            using Oceananigans.Advection: stencil_sum, $stencil, $weno_weights, $biased_p
            import Oceananigans.Advection: $interpolate_func

            @inline function $interpolate_func(i, j, k, grid, 
                                                scheme::WENOVectorInvariant{N, FT, XT, YT, ZT}, 
                                                ψ, idx, loc, VI::Type{VelocityStencil}, args...) where {N, FT, XT, YT, ZT}

                @inbounds begin
                    ψₜ = $stencil(i, j, k, scheme, ψ, grid, args...)
                    w = $weno_weights((i, j, k, idx), scheme, Val($val), VI, grid, args...)
                end
                return stencil_sum(scheme, ψₜ, w, $biased_p, $cT, $val, idx, loc)
            end
        end
    end
end

using Oceananigans.Operators
using Oceananigans.Operators: ℑxyᶠᶠᵃ, div_xyᶜᶜᶜ
using Oceananigans.Advection: ZWENO, Cl, Cr

import Oceananigans.Advection: tangential_left_stencil_u, tangential_left_stencil_v

@inline interpolation_u(i, idx) = ifelse(i == idx, ℑxᶜᵃᵃ, ℑyᵃᶠᵃ)
@inline interpolation_v(j, jdx) = ifelse(j == jdx, ℑyᵃᶜᵃ, ℑxᶠᵃᵃ)

@inline tangential_left_stencil_u(i, j, k, idx, scheme, ::Val{1}, u) = @inbounds left_stencil_x(i, j, k, scheme, interpolation_u(i, idx), u)
@inline tangential_left_stencil_u(i, j, k, idx, scheme, ::Val{2}, u) = @inbounds left_stencil_y(i, j, k, scheme, interpolation_u(j, idx), u)
@inline tangential_left_stencil_v(i, j, k, idx, scheme, ::Val{1}, v) = @inbounds left_stencil_x(i, j, k, scheme, interpolation_v(i, idx), v)
@inline tangential_left_stencil_v(i, j, k, idx, scheme, ::Val{2}, v) = @inbounds left_stencil_y(i, j, k, scheme, interpolation_v(j, idx), v)

@inline tangential_right_stencil_u(i, j, k, idx, scheme, ::Val{1}, u) = @inbounds right_stencil_x(i, j, k, scheme, interpolation_u(i, idx), u)
@inline tangential_right_stencil_u(i, j, k, idx, scheme, ::Val{2}, u) = @inbounds right_stencil_y(i, j, k, scheme, interpolation_u(j, idx), u)
@inline tangential_right_stencil_v(i, j, k, idx, scheme, ::Val{1}, v) = @inbounds right_stencil_x(i, j, k, scheme, interpolation_v(i, idx), v)
@inline tangential_right_stencil_v(i, j, k, idx, scheme, ::Val{2}, v) = @inbounds right_stencil_y(i, j, k, scheme, interpolation_v(j, idx), v)

# Calculating Dynamic WENO Weights (wᵣ), either with JS weno, Z weno or VectorInvariant WENO
for (side, coeff) in zip([:left, :right], (:Cl, :Cr))
    biased_weno_weights = Symbol(side, :_biased_weno_weights)
    biased_β = Symbol(side, :_biased_β)
    
    divergence_stencil = Symbol(:divergence_, side, :_stencil)
    
    tangential_stencil_u = Symbol(:tangential_, side, :_stencil_u)
    tangential_stencil_v = Symbol(:tangential_, side, :_stencil_v)

    @eval begin
        using Oceananigans.Advection: beta_loop, global_smoothness_indicator, zweno_alpha_loop, js_alpha_loop
        using Oceananigans.Advection: $biased_β

        import Oceananigans.Advection:  $biased_weno_weights

        @inline function $biased_weno_weights(ijk, scheme::WENO{N, FT}, dir, ::Type{VelocityStencil}, u, v) where {N, FT}
            @inbounds begin
                i, j, k, idx = ijk
            
                uₛ = $tangential_stencil_u(i, j, k, idx, scheme, dir, u)
                vₛ = $tangential_stencil_v(i, j, k, idx, scheme, dir, v)
            
                βᵤ = beta_loop(scheme, uₛ, $biased_β)
                βᵥ = beta_loop(scheme, vₛ, $biased_β)

                β  = beta_sum(scheme, βᵤ, βᵥ)

                if scheme isa ZWENO
                    τ = global_smoothness_indicator(Val(N), β)
                    α = zweno_alpha_loop(scheme, β, τ, $coeff, FT)
                else
                    α = js_alpha_loop(scheme, β, $coeff, FT)
                end
                return α ./ sum(α)
            end
        end
    end
end
