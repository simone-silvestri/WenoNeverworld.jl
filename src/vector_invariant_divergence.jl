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

import Oceananigans.Advection: vertical_vorticity_U, vertical_vorticity_V

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

    @inbounds û = u[i, j, k]
    ∂uᴸ =  _left_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, ∂xᶜᶜᶜ, u)
    ∂uᴿ = _right_biased_interpolate_xᶠᵃᵃ(i, j, k, grid, scheme, ∂xᶜᶜᶜ, u)

    @inbounds v̂ = ℑxᶠᵃᵃ(i, j, k, grid, ℑyᵃᶜᵃ, Δx_qᶜᶠᶜ, v) / Δxᶠᶜᶜ(i, j, k, grid) 
    ∂vᴸ =  _left_biased_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme, ∂xᶠᶠᶜ, v)
    ∂vᴿ = _right_biased_interpolate_yᵃᶜᵃ(i, j, k, grid, scheme, ∂xᶠᶠᶜ, v)

    return upwind_biased_product(û, ∂uᴸ, ∂uᴿ) + upwind_biased_product(û, ∂vᴸ, ∂vᴿ)
end

@inline function bernoulli_head_V(i, j, k, grid, scheme::WENOVectorInvariant, u, v)

    @inbounds û = ℑyᵃᶠᵃ(i, j, k, grid, ℑxᶜᵃᵃ, Δy_qᶠᶜᶜ, u) / Δyᶜᶠᶜ(i, j, k, grid)
    ∂uᴸ =  _left_biased_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme, ∂yᶠᶠᶜ, u)
    ∂uᴿ = _right_biased_interpolate_xᶜᵃᵃ(i, j, k, grid, scheme, ∂yᶠᶠᶜ, u)

    @inbounds v̂ = v[i, j, k]
    ∂vᴸ =  _left_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, ∂yᶜᶜᶜ, v)
    ∂vᴿ = _right_biased_interpolate_yᵃᶠᵃ(i, j, k, grid, scheme, ∂yᶜᶜᶜ, v)

    return upwind_biased_product(û, ∂uᴸ, ∂uᴿ) + upwind_biased_product(û, ∂vᴸ, ∂vᴿ)
end


# @inline function vertical_vorticity_U(i, j, k, grid, scheme::WENOVectorInvariant{N, FT, XT, YT, ZT, VI}, u, v) where {N, FT, XT, YT, ZT, VI}
#     v̂  =  ℑxᶠᵃᵃ(i, j, k, grid, ℑyᵃᶜᵃ, Δx_qᶜᶠᶜ, v) / Δxᶠᶜᶜ(i, j, k, grid) 
#     ζᴸ =  _multi_dimensional_reconstruction_x(i, j, k, grid, scheme,  _left_biased_interpolate_yᵃᶜᵃ, ζ₃ᶠᶠᶜ, VI, u, v)
#     ζᴿ =  _multi_dimensional_reconstruction_x(i, j, k, grid, scheme, _right_biased_interpolate_yᵃᶜᵃ, ζ₃ᶠᶠᶜ, VI, u, v)
#     return - upwind_biased_product(v̂, ζᴸ, ζᴿ) 
# end

# @inline function vertical_vorticity_V(i, j, k, grid, scheme::WENOVectorInvariant{N, FT, XT, YT, ZT, VI}, u, v) where {N, FT, XT, YT, ZT, VI}
#     û  =  ℑyᵃᶠᵃ(i, j, k, grid, ℑxᶜᵃᵃ, Δy_qᶠᶜᶜ, u) / Δyᶜᶠᶜ(i, j, k, grid)
#     ζᴸ =  _multi_dimensional_reconstruction_y(i, j, k, grid, scheme,  _left_biased_interpolate_xᶜᵃᵃ, ζ₃ᶠᶠᶜ, VI, u, v)
#     ζᴿ =  _multi_dimensional_reconstruction_y(i, j, k, grid, scheme, _right_biased_interpolate_xᶜᵃᵃ, ζ₃ᶠᶠᶜ, VI, u, v)
#     return + upwind_biased_product(û, ζᴸ, ζᴿ) 
# end

# @inline _multi_dimensional_reconstruction_x(args...) = multi_dimensional_reconstruction_x(args...)
# @inline _multi_dimensional_reconstruction_y(args...) = multi_dimensional_reconstruction_y(args...)


# const ε = 1e-6
# const two_32 = Int32(2)

# const σ⁺ = 214/80
# const σ⁻ =  67/40

# ## Figure them out!
# const γ₀¹  =   9.0 / 80 / σ⁺
# const γ₁¹  =  49.0 / 20 / σ⁺
# const γ₂¹  =   9.0 / 80 / σ⁺

# const γ₀²⁺ =   9.0 / 80 / σ⁺
# const γ₁²⁺ =  49.0 / 20 / σ⁺
# const γ₂²⁺ =   9.0 / 80 / σ⁺

# const γ₀²⁻ =   9.0 / 40 / σ⁻
# const γ₁²⁻ =  49.0 / 40 / σ⁻
# const γ₂²⁻ =   9.0 / 40 / σ⁻

# const γ₀³  =   9.0 / 40 / σ⁻
# const γ₁³  =  49.0 / 40 / σ⁻
# const γ₂³  =   9.0 / 40 / σ⁻

# ## Figure them out!
# const a₀¹ = (2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) / 60
# const a₁¹ = (2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) / 60
# const a₂¹ = (2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) / 60

# const a₀² = (2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) / 60
# const a₁² = (2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) / 60
# const a₂² = (2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) / 60

# const a₀³ = (2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) / 60
# const a₁³ = (2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) / 60
# const a₂³ = (2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) / 60

# @inline function multi_dimensional_reconstruction_y(i, j, k, grid, scheme, _interpolate_x, f::Function, VI, args...)

#     Q₋₂ = _interpolate_x(i, j, k, grid, scheme, f, VI, args...)
#     Q₋₁ = _interpolate_x(i, j, k, grid, scheme, f, VI, args...)
#     Q₀  = _interpolate_x(i, j, k, grid, scheme, f, VI, args...)
#     Q₊₁ = _interpolate_x(i, j, k, grid, scheme, f, VI, args...)
#     Q₊₂ = _interpolate_x(i, j, k, grid, scheme, f, VI, args...)

#     S₀ = (Q₋₂, Q₋₁, Q₀)
#     S₁ = (Q₋₁, Q₀ , Q₊₁)
#     S₂ = (Q₀ , Q₊₁, Q₊₂)

#     q = zeros(3)

#     @unroll for j in 1:3


#     end

# end