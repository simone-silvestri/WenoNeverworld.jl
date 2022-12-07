
@inline _multi_dimensional_reconstruction_x(args...) = multi_dimensional_reconstruction_x(args...)
@inline _multi_dimensional_reconstruction_y(args...) = multi_dimensional_reconstruction_y(args...)

const ε = 1e-6
const two_32 = Int32(2)

## Figure them out!
const γ₀¹  = (1008 + 71 * sqrt(15)) / 5240
const γ₁¹  =  408 / 655
const γ₂¹  = (1008 - 71 * sqrt(15)) / 5240

const σ⁺ = 214/80
const σ⁻ =  67/40

const γ₀²⁺ =   9.0 / 80 / σ⁺
const γ₁²⁺ =  49.0 / 20 / σ⁺
const γ₂²⁺ =   9.0 / 80 / σ⁺

const γ₀²⁻ =   9.0 / 40 / σ⁻
const γ₁²⁻ =  49.0 / 40 / σ⁻
const γ₂²⁻ =   9.0 / 40 / σ⁻

const γ₀³  = (1008 - 71 * sqrt(15)) / 5240
const γ₁³  =  408 / 655
const γ₂³  = (1008 + 71 * sqrt(15)) / 5240

## Figure them out!
const a₀¹ = ( 2 - 3*sqrt(15), -4 + 12*sqrt(15), 62 - 9 * sqrt(15)) ./ 60
const a₁¹ = ( 2 + 3*sqrt(15),               56,  2 - 3 * sqrt(15)) ./ 60
const a₂¹ = (62 + 9*sqrt(15), -4 - 12*sqrt(15),  2 + 3 * sqrt(15)) ./ 60

const a₀² = (-1, 2,  23) ./ 24
const a₁² = (-1, 26, -1) ./ 24
const a₂² = (23, 2,  -1) ./ 24

const a₀³ = ( 2 + 3*sqrt(15), -4 - 12*sqrt(15), 62 + 9 * sqrt(15)) ./ 60
const a₁³ = ( 2 - 3*sqrt(15),               56,  2 + 3 * sqrt(15)) ./ 60
const a₂³ = (62 - 9*sqrt(15), -4 + 12*sqrt(15),  2 - 3 * sqrt(15)) ./ 60

const a₀ˢ = (a₀¹, a₀², a₀³)
const a₁ˢ = (a₁¹, a₁², a₁³)
const a₂ˢ = (a₂¹, a₂², a₂³)

const aₜˢ = (a₀ˢ, a₁ˢ, a₂ˢ)

@inline function multi_dimensional_reconstruction_x(i, j, k, grid, scheme, _interpolate_y, f::Function, VI, args...)

    FT = eltype(grid)

    Q₋₂ = _interpolate_y(i-2, j, k, grid, scheme, f, VI, args...)
    Q₋₁ = _interpolate_y(i-1, j, k, grid, scheme, f, VI, args...)
    Q₀  = _interpolate_y(i,   j, k, grid, scheme, f, VI, args...)
    Q₊₁ = _interpolate_y(i+1, j, k, grid, scheme, f, VI, args...)
    Q₊₂ = _interpolate_y(i+2, j, k, grid, scheme, f, VI, args...)

    S₀ = (Q₋₂, Q₋₁, Q₀)
    S₁ = (Q₋₁, Q₀ , Q₊₁)
    S₂ = (Q₀ , Q₊₁, Q₊₂)

    return fifth_order_weno_reconstruction(FT, S₀, S₁, S₂)
end

@inline function multi_dimensional_reconstruction_y(i, j, k, grid, scheme, _interpolate_x, f::Function, VI, args...)

    FT = eltype(grid)

    Q₋₂ = _interpolate_x(i, j-2, k, grid, scheme, f, VI, args...)
    Q₋₁ = _interpolate_x(i, j-1, k, grid, scheme, f, VI, args...)
    Q₀  = _interpolate_x(i, j,   k, grid, scheme, f, VI, args...)
    Q₊₁ = _interpolate_x(i, j+1, k, grid, scheme, f, VI, args...)
    Q₊₂ = _interpolate_x(i, j+2, k, grid, scheme, f, VI, args...)

    S₀ = (Q₋₂, Q₋₁, Q₀)
    S₁ = (Q₋₁, Q₀ , Q₊₁)
    S₂ = (Q₀ , Q₊₁, Q₊₂)

    return fifth_order_weno_reconstruction(FT, S₀, S₁, S₂)
end

function fifth_order_weno_reconstruction(FT, S₀, S₁, S₂)

    q̂₀¹ = sum(a₀¹ .* S₀)
    q̂₁¹ = sum(a₁¹ .* S₁)
    q̂₂¹ = sum(a₂¹ .* S₂)

    q̂₀² = sum(a₀² .* S₀)
    q̂₁² = sum(a₁² .* S₁)
    q̂₂² = sum(a₂² .* S₂)
    
    q̂₀³ = sum(a₀³ .* S₀)
    q̂₁³ = sum(a₁³ .* S₁)
    q̂₂³ = sum(a₂³ .* S₂)

    β₀ =   left_biased_β_constant(FT, S₀)
    β₁ = center_biased_β_constant(FT, S₁)
    β₂ =  right_biased_β_constant(FT, S₂)

    w₀¹, w₁¹, w₂¹ = centered_reconstruction_weights¹(FT, β₀, β₁, β₂)
    w₀³, w₁³, w₂³ = centered_reconstruction_weights³(FT, β₀, β₁, β₂)

    w₀²⁺, w₁²⁺, w₂²⁺, w₀²⁻, w₁²⁻, w₂²⁻ = centered_reconstruction_weights²(FT, β₀, β₁, β₂)

    q¹ = w₀¹ * q̂₀¹ + w₁¹ * q̂₁¹ + w₂¹ * q̂₂¹
    q³ = w₀³ * q̂₀³ + w₁³ * q̂₁³ + w₂³ * q̂₂³

    q²⁺ = w₀²⁺ * q̂₀² + w₁²⁺ * q̂₁² + w₂²⁺ * q̂₂²
    q²⁻ = w₀²⁻ * q̂₀² + w₁²⁻ * q̂₁² + w₂²⁻ * q̂₂²

    q² = σ⁺ * q²⁺ - σ⁻ * q²⁻

    return sqrt(π) * (q¹ / 6 + 2 * q² / 3 + q³ / 6)
end

@inline   left_biased_β_constant(FT, ψ) = @inbounds FT(13/12) * (ψ[1] - 2ψ[2] + ψ[3])^two_32 + FT(1/4) * ( ψ[1] - 4ψ[2] + 3ψ[3])^two_32
@inline center_biased_β_constant(FT, ψ) = @inbounds FT(13/12) * (ψ[1] - 2ψ[2] + ψ[3])^two_32 + FT(1/4) * ( ψ[1]         -  ψ[3])^two_32
@inline  right_biased_β_constant(FT, ψ) = @inbounds FT(13/12) * (ψ[1] - 2ψ[2] + ψ[3])^two_32 + FT(1/4) * (3ψ[1] - 4ψ[2] +  ψ[3])^two_32

@inline function centered_reconstruction_weights²(FT, β₀, β₁, β₂)

    α₀⁺ = FT(γ₀²⁺) / (β₀ + FT(ε))^two_32
    α₁⁺ = FT(γ₁²⁺) / (β₁ + FT(ε))^two_32
    α₂⁺ = FT(γ₂²⁺) / (β₂ + FT(ε))^two_32

    α₀⁻ = FT(γ₀²⁻) / (β₀ + FT(ε))^two_32
    α₁⁻ = FT(γ₁²⁻) / (β₁ + FT(ε))^two_32
    α₂⁻ = FT(γ₂²⁻) / (β₂ + FT(ε))^two_32

    Σα⁺ = α₀⁺ + α₁⁺ + α₂⁺
    w₀⁺ = α₀⁺ / Σα⁺
    w₁⁺ = α₁⁺ / Σα⁺
    w₂⁺ = α₂⁺ / Σα⁺

    Σα⁻ = α₀⁻ + α₁⁻ + α₂⁻
    w₀⁻ = α₀⁻ / Σα⁻
    w₁⁻ = α₁⁻ / Σα⁻
    w₂⁻ = α₂⁻ / Σα⁻

    return w₀⁺, w₁⁺, w₂⁺, w₀⁻, w₁⁻, w₂⁻
end

@inline function centered_reconstruction_weights¹(FT, β₀, β₁, β₂)

    α₀ = FT(γ₀¹) / (β₀ + FT(ε))^two_32
    α₁ = FT(γ₁¹) / (β₁ + FT(ε))^two_32
    α₂ = FT(γ₂¹) / (β₂ + FT(ε))^two_32

    Σα = α₀ + α₁ + α₂
    w₀ = α₀ / Σα
    w₁ = α₁ / Σα
    w₂ = α₂ / Σα
    
    return w₀, w₁, w₂
end

@inline function centered_reconstruction_weights³(FT, β₀, β₁, β₂)

    α₀ = FT(γ₀³) / (β₀ + FT(ε))^two_32
    α₁ = FT(γ₁³) / (β₁ + FT(ε))^two_32
    α₂ = FT(γ₂³) / (β₂ + FT(ε))^two_32

    Σα = α₀ + α₁ + α₂
    w₀ = α₀ / Σα
    w₁ = α₁ / Σα
    w₂ = α₂ / Σα

    return w₀, w₁, w₂
end