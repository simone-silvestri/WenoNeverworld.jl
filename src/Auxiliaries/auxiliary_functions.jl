using WenoNeverworld.Constants

""" utility profiles (atan, exponential, and parabolic) """
@inline exponential_profile(z; Δ = Constants.ΔB, Lz = Constants.Lz, h = Constants.h) = ( Δ * (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) )
@inline parabolic_scaling(y) = - 1 / Constants.max_latitude^2 * y^2 + 1
@inline initial_buoyancy_parabola(x, y, z) = exponential_profile(z) * parabolic_scaling(y) 

""" 
    function cubic_interpolate(x, x1, x2, y1, y2, d1, d2)

returns a cubic function between points `(x1, y1)` and `(x2, y2)` with derivative `d1` and `d2`
"""
@inline function cubic_interpolate(x; x₁, x₂, y₁, y₂, d₁ = 0, d₂ = 0)
    A = [ x₁^3 x₁^2 x₁ 1.0
          x₂^3 x₂^2 x₂ 1.0
          3*x₁^2 2*x₁ 1.0 0.0
          3*x₂^2 2*x₂ 1.0 0.0]
          
    b = [y₁, y₂, d₁, d₂]

    coeff = A \ b

    return coeff[1] * x^3 + coeff[2] * x^2 + coeff[3] * x + coeff[4]
end
