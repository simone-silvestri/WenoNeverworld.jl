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

using Oceananigans.Operators: ℑxyᶠᶠᵃ, div_xyᶜᶜᶜ

@inline divergence_left_stencil(i, j, k, grid, scheme, ::Val{1}, u, v) = @inbounds left_stencil_x(i, j, k, scheme, ℑxyᶠᶠᵃ, grid, div_xyᶜᶜᶜ, u, v)
@inline divergence_left_stencil(i, j, k, grid, scheme, ::Val{2}, u, v) = @inbounds left_stencil_y(i, j, k, scheme, ℑxyᶠᶠᵃ, grid, div_xyᶜᶜᶜ, u, v)

@inline divergence_right_stencil(i, j, k, grid, scheme, ::Val{1}, u, v) = @inbounds right_stencil_x(i, j, k, scheme, ℑxyᶠᶠᵃ, grid, div_xyᶜᶜᶜ, u, v)
@inline divergence_right_stencil(i, j, k, grid, scheme, ::Val{2}, u, v) = @inbounds right_stencil_y(i, j, k, scheme, ℑxyᶠᶠᵃ, grid, div_xyᶜᶜᶜ, u, v)

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
