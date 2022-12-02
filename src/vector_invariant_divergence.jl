import Oceananigans.Advection: vertical_vorticity_U, vertical_vorticity_V

using Oceananigans.Operators
using Oceananigans.Advection: WENOVectorInvariant, 
		    _left_biased_interpolate_yᵃᶜᵃ,
		    _left_biased_interpolate_xᶜᵃᵃ,
		   _right_biased_interpolate_yᵃᶜᵃ,
		   _right_biased_interpolate_xᶜᵃᵃ,
		    	    upwind_biased_product

using Oceananigans.Advection: SmoothnessStencil

struct DivergenceStencil <:SmoothnessStencil end

const WENOVectorInvariantDiv{N, FT, XT, YT, ZT, VI, WF, PP}  = 
      WENO{N, FT, XT, YT, ZT, VI, WF, PP} where {N, FT, XT, YT, ZT, VI<:DivergenceStencil, WF, PP}

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
                                               ψ, idx, loc, VI::Type{DivergenceStencil}, args...) where {N, FT, XT, YT, ZT}

                @inbounds begin
                    ψₜ = $stencil(i, j, k, scheme, ψ, grid, args...)
                    w = $weno_weights((i, j, k), scheme, Val($val), VI, grid, args...)
                    return stencil_sum(scheme, ψₜ, w, $biased_p, $cT, $val, idx, loc)
                end
            end
        end
    end
end

using Oceananigans.Operators

@inline div_xyᶠᶠᶜ(i, j, k, grid, u, v) = 
    1 / Azᶜᶜᶜ(i, j, k, grid) * (δxᶠᵃᵃ(i, j, k, grid, Δy_qᶜᶠᶜ, ℑxyᶜᶠᵃ, u) +
                                δyᵃᶠᵃ(i, j, k, grid, Δx_qᶠᶜᶜ, ℑxyᶠᶜᵃ, v))

using Oceananigans.Operators: ℑxyᶠᶠᵃ, div_xyᶜᶜᶜ

@inline divergence_left_stencil(i, j, k, grid, scheme, ::Val{1}, u, v) = @inbounds left_stencil_x(i, j, k, scheme, div_xyᶜᶜᶜ, grid, u, v)
@inline divergence_left_stencil(i, j, k, grid, scheme, ::Val{2}, u, v) = @inbounds left_stencil_y(i, j, k, scheme, div_xyᶜᶜᶜ, grid, u, v)

@inline divergence_right_stencil(i, j, k, grid, scheme, ::Val{1}, u, v) = @inbounds right_stencil_x(i, j, k, scheme, div_xyᶜᶜᶜ, grid, u, v)
@inline divergence_right_stencil(i, j, k, grid, scheme, ::Val{2}, u, v) = @inbounds right_stencil_y(i, j, k, scheme, div_xyᶜᶜᶜ, grid, u, v)

using Oceananigans.Advection: ZWENO, Cl, Cr

# Calculating Dynamic WENO Weights (wᵣ), either with JS weno, Z weno or VectorInvariant WENO
for (side, coeff) in zip([:left, :right], (:Cl, :Cr))
    biased_weno_weights = Symbol(side, :_biased_weno_weights)
    biased_β = Symbol(side, :_biased_β)
    
    divergence_stencil = Symbol(:divergence_, side, :_stencil)
    
    @eval begin

        using Oceananigans.Advection: beta_loop, global_smoothness_indicator, zweno_alpha_loop, js_alpha_loop
        using Oceananigans.Advection: $biased_β

        import Oceananigans.Advection:  $biased_weno_weights

        @inline function $biased_weno_weights(ijk, scheme::WENO{N, FT}, dir, ::Type{DivergenceStencil}, grid, u, v) where {N, FT}
            @inbounds begin
                i, j, k = ijk
            
                δ = $divergence_stencil(i, j, k, grid, scheme, dir, u, v)
                β = beta_loop(scheme, δ, $biased_β)

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

using Oceananigans.Advection: VorticityStencil

for bias in (:left_biased, :right_biased)
    for (d, dir) in zip((:x, :y), (:xᶜᵃᵃ, :yᵃᶜᵃ))
        interp     = Symbol(bias, :_interpolate_, dir)
        alt_interp = Symbol(:_, interp)

        near_horizontal_boundary = Symbol(:near_, d, :_horizontal_boundary_, bias)

        @eval begin
            using Oceananigans.ImmersedBoundaries: $near_horizontal_boundary
            using Oceananigans.Advection: $interp
            import Oceananigans.Advection: $alt_interp

            # Conditional Interpolation for DivergenceStencil WENO vector invariant scheme
            @inline $alt_interp(i, j, k, ibg::ImmersedBoundaryGrid, scheme::WENOVectorInvariantDiv, ζ, ::Type{DivergenceStencil}, u, v) =
                ifelse($near_horizontal_boundary(i, j, k, ibg, scheme),
                    $alt_interp(i, j, k, ibg, scheme, ζ, VorticityStencil, u, v),
                    $interp(i, j, k, ibg, scheme, ζ, DivergenceStencil, u, v))
        end
    end
end

import Oceananigans.Advection: WENO
using Oceananigans.Advection: compute_reconstruction_coefficients

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
        buffer_scheme   = WENO(vector_invariant, FT; grid, order = order - 2, zweno, bounds)
        advecting_velocity_scheme = Centered(FT; grid, order = order - 1)
    end

    return WENO{N, FT, VI, zweno}(weno_coefficients..., bounds, buffer_scheme, advecting_velocity_scheme)
end

import Oceananigans.Advection: vertical_advection_U, vertical_advection_V
import Oceananigans.Advection: bernoulli_head_U, bernoulli_head_V
import Oceananigans.Advection: U_dot_∇u, U_dot_∇v

using Oceananigans.Advection: vertical_vorticity_U, vertical_vorticity_V

@inline U_dot_∇u(i, j, k, grid, scheme::WENOVectorInvariant, U) = (
    + vertical_vorticity_U(i, j, k, grid, scheme, U.u, U.v)  # Vertical relative vorticity term
    + vertical_advection_U(i, j, k, grid, scheme, U)  # Horizontal vorticity / vertical advection term
    + bernoulli_head_U(i, j, k, grid, scheme, U.u, U.v))     # Bernoulli head term
    
@inline U_dot_∇v(i, j, k, grid, scheme::WENOVectorInvariant, U) = (
    + vertical_vorticity_V(i, j, k, grid, scheme, U.u, U.v)  # Vertical relative vorticity term
    + vertical_advection_V(i, j, k, grid, scheme, U)  # Horizontal vorticity / vertical advection term
    + bernoulli_head_V(i, j, k, grid, scheme, U.u, U.v))     # Bernoulli head term

using Oceananigans.Advection: _advective_momentum_flux_Wu, 
                              _advective_momentum_flux_Wv,
                              _advective_momentum_flux_Uu,
                              _advective_momentum_flux_Vv

using Oceananigans.Advection:  _left_biased_interpolate_xᶠᵃᵃ,
                               _left_biased_interpolate_yᵃᶠᵃ,
                              _right_biased_interpolate_xᶠᵃᵃ,
                              _right_biased_interpolate_yᵃᶠᵃ

@inline function vertical_advection_U(i, j, k, grid, scheme::WENOVectorInvariant, U)

    u, v, W = U

    @inbounds û = u[i, j, k] 
    δᴸ    =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, div_xyᶜᶜᶜ, VorticityStencil, u, v)
    δᴿ    = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, div_xyᶜᶜᶜ, VorticityStencil, u, v)
    δflux = upwind_biased_product(û, δᴸ, δᴿ) 
    wterm =  1/Vᶠᶜᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wu, scheme, W, u)

    return wterm + δflux
end

@inline function vertical_advection_V(i, j, k, grid, scheme::WENOVectorInvariant, U)

    u, v, W = U

    @inbounds v̂ = v[i, j, k] 
    δᴸ    =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, div_xyᶜᶜᶜ, VorticityStencil, u, v)
    δᴿ    = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, div_xyᶜᶜᶜ, VorticityStencil, u, v)
    δflux = upwind_biased_product(v̂, δᴸ, δᴿ) 
    wterm =  1/Vᶜᶠᶜ(i, j, k, grid) * δzᵃᵃᶜ(i, j, k, grid, _advective_momentum_flux_Wv, scheme, W, v)

    return wterm + δflux
end

@inline function bernoulli_head_U(i, j, k, grid, scheme::WENOVectorInvariant, u, v)

    uterm =  1/Vᶠᶜᶜ(i, j, k, grid) * δxᶠᵃᵃ(i, j, k, grid, _advective_momentum_flux_Uu, scheme, u, u)
    vterm = ∂yᶜᶠᶜ(i, j, k, grid, ℑyᵃᶜᵃ, ϕ², v)

    return 0.5 * (uterm + vterm)
end

@inline function bernoulli_head_V(i, j, k, grid, scheme::WENOVectorInvariant, u, v)

    vterm =  1/Vᶜᶠᶜ(i, j, k, grid) * δyᵃᶠᵃ(i, j, k, grid, _advective_momentum_flux_Vv, scheme, v, v)
    uterm = ∂xᶠᶜᶜ(i, j, k, grid, ℑxᶜᵃᵃ, ϕ², u)

    return 0.5 * (uterm + vterm)
end