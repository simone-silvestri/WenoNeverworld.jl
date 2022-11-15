function cubic_profile(x, x1, x2, y1, y2, d1, d2)
    A = [ x1^3 x1^2 x1 1.0
          x2^3 x2^2 x2 1.0
          3*x1^2 2*x1 1.0 0.0
          3*x2^2 2*x2 1.0 0.0]
          
    b = [y1, y2, d1, d2]

    coeff = A \ b

    return coeff[1] * x^3 + coeff[2] * x^2 + coeff[3] * x + coeff[4]
end

## define the coasts
function coastal_ridge_x(x) 
    if x < 0.5
        return 0.0
    elseif x < 2.5 
        return -200.0
    elseif x < 5.0
        return cubic_profile(x, 2.5, 5.0, -200.0, -4000, 0.0, 0.0)
    else        
        return -4000.0 
    end
end

function coastal_ridge_y(x) 
    if x < 0.5
        return cubic_profile(x, 0.0, 0.5, 0.0, -200, 0.0, 0.0)
    elseif x < 2.5 
        return -200.0
    elseif x < 5.0
        return cubic_profile(x, 2.5, 5.0, -200.0, -4000, 0.0, 0.0)
    else        
        return -4000.0 
    end
end

function bottom_ridge_x(x)
    if x < 20
        return -4000
    elseif x < 29
        return cubic_profile(x, 20.0, 29.0, -4000, -2000, 0.0, 0.0)
    elseif x < 31
        return -2000.0
    else 
        return -4000.0
    end
end

function bottom_ridge_xy(x, y)
    if y > - 30
        return bottom_ridge_x(x)
    elseif y > -50
        return cubic_profile(y, -30, -50, bottom_ridge_x(x), -4000, 0.0, 0.0)
    else
        return -4000.0
    end
end
        
function scotia_arc(x, y)
    radius = sqrt(x^2 + (y + 50)^2)
    if radius < 8
        return -4000.0
    elseif radius < 9
        return cubic_profile(radius, 8.0, 9.0, -4000.0, -2000.0, 0.0, 0.0)
    elseif radius < 11
        return -2000.0
    elseif radius < 12
        return cubic_profile(radius, 11.0, 12.0, -2000.0, -4000.0, 0.0, 0.0)
    else
        return -4000.0
    end
end

function limit_last_cells(y)
    if y < -0.5 || y > -69.5
        return -4000.0
    else
        return 0.0
    end
end

function bathymetry(x, y) 
    if x < 5 || x > 55
        if x < 0 
           x = 0.0
        end
        if x > 60
           x = 60.0
        end
        if y > -59 && y < -41 
            return  max(scotia_arc(x, y), 
                       coastal_ridge_x(sqrt(x^2 + (y + 59)^2)),
                       coastal_ridge_x(sqrt(x^2 + (y + 41)^2)), 
                       coastal_ridge_x(sqrt((60 - x)^2 + (y + 59)^2)),
                       coastal_ridge_x(sqrt((60 - x)^2 + (y + 41)^2)))
        else
            return max(coastal_ridge_x(x), 
                       coastal_ridge_x(60 - x), 
                       bottom_ridge_xy(x, y), 
                       bottom_ridge_xy(60 - x, y), 
                       scotia_arc(x, y))
        end
    else
        return max(coastal_ridge_x(x), 
                   coastal_ridge_x(60 - x), 
                   bottom_ridge_xy(x, y), 
                   bottom_ridge_xy(60 - x, y), 
                   scotia_arc(x, y))
    end
end

@inline function wind_stress(y, mid_wind)
    if y < -45
        return cubic_profile(y, -70.0, -45.0, 0.0, 0.2, 0.0, 0.0)
    elseif y < -15
        return cubic_profile(y, -45.0, -15.0, 0.2, -0.1, 0.0, 0.0)
    elseif y < 0
        return cubic_profile(y, -15.0, 0.0, -0.1, mid_wind, 0.0, 0.0)
    elseif y < 15
        return cubic_profile(y, 0.0, 15.0, mid_wind, -0.1, 0.0, 0.0)
    elseif y < 45
        return cubic_profile(y, 15.0, 45.0, -0.1, 0.1, 0.0, 0.0)
    else
        return cubic_profile(y, 45.0, 70.0, 0.1, 0.0, 0.0, 0.0)
    end
end

const Lz   = 4000.0
const Ly   = 70
const h    = 1000.0
const ΔB   = 6.0e-2 
const fact = 5.0

@inline atan_scaling(y) = (atan(fact*((Ly + y)/Ly - 0.5)) / atan(fact * 0.5) + 1) /2
@inline initial_buoyancy(x, y, z) = ( ΔB * (exp(z / h) - exp( - Lz / h)) / (1 - exp( - Lz / h)) ) * atan_scaling(y)
