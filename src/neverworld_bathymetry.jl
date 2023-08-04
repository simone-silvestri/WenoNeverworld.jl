# The bathymetry is defined for a latitude range of -70 ≤ φ ≤ 70
# and a longitude range of 0 ≤ λ ≤ 60

Base.@kwdef struct ShelfParameters
    coast_length::Float64 = 0.5
    side_length::Float64  = 2.5
    size::Float64         = 2.5
    depth::Float64        = 200
end

Base.@kwdef struct RidgeParameters
    side_length::Float64 = 9
    top_length::Float64  = 2
    mid_point::Float64   = 30
    depth::Float64       = 2000
end

Base.@kwdef struct ScotiaArcParameters
    left_outer_radius::Float64  = 8
    left_inner_radius::Float64  = 9
    right_inner_radius::Float64 = 11
    right_outer_radius::Float64 = 12
    mid_point::Float64          = 50
    depth::Float64              = 2000
end

Base.@kwdef struct NeverWorldBathymetryParameters
    shelves::ShelfParameters        = ShelfParameters()
    scotia_arc::ScotiaArcParameters = ScotiaArcParameters()
    channel_south_edge::Float64     = - 59
    channel_north_edge::Float64     = - 41
    bottom::Float64                 = - 4000
    ridge::Nothing                  = nothing
end
    
## define the coasts
function coastal_shelf_x(x, params, bottom) 

    coast       = params.coast_length
    length      = params.length
    side_length = params.side_length + length
    depth       = - params.shelf_depth

    if x < coast
        return 0.0
    elseif x < length
        return depth
    elseif x < side_length
        return cubic_profile(x, x₁ = length, x₂ = side_length, y₁ = depth, y₂ = bottom)
    else        
        return bottom
    end
end

function sharp_coast_x(x, params, bottom) 

    coast  = params.coast_length
    bottom = - params.bottom_depth

    if x < coast
        return 0.0
    else        
        return bottom
    end
end

function coastal_shelf_y(x, params, bottom) 

    coast       = params.coast_length
    length      = params.length
    side_length = params.side_length + length
    depth       = - params.shelf_depth

    if x < coast
        return cubic_profile(x, x₁ = 0.0, x₂ = coast, y₁ = 0.0, y₂ = depth)
    elseif x < length
        return depth
    elseif x < side_length
        return cubic_profile(x, x₁ = length, x₂ = side_length, y₁ = depth, y₂ = bottom)
    else        
        return bottom
    end
end

bottom_ridge_x(x, ::Nothing, bottom) = bottom

# Bottom ridge
function bottom_ridge_x(x, params, bottom)

    center     = params.mid_point_longitude
    
    top_left   = center - params.top_length/2
    top_right  = center + params.top_length/2

    bot_left  = center - params.top_length/2 - params.ridge_side_length
    bot_right = center + params.top_length/2 + params.ridge_side_length

    depth = - params.ridge_depth

    if x < bot_left
        return bottom
    elseif x < top_left
        return cubic_profile(x, x₁ = bot_left, x₂ = top_left, y₁ = bottom, y₂ = depth)
    elseif x < top_right
        return depth
    elseif x < bot_right
        return cubic_profile(x, x₁ = top_right, x₂ = bot_right, y₁ = depth, y₂ = ridge)
    else 
        return bottom
    end
end
        
# Scotia arc
function scotia_arc(x, y, params, bottom)
    
    left_inner_radius = params.left_inner_radius
    left_outer_radius = params.left_outer_radius
    right_inner_radius = params.right_inner_radius
    right_outer_radius = params.right_outer_radius
    mid_point = params.mid_point
    depth = - params.depth

    radius = sqrt(x^2 + (y + mid_point)^2)
    if radius < left_inner_radius
        return bottom
    elseif radius < left_outer_radius
        return cubic_profile(radius, x₁ = left_outer_radius, x₂ = left_inner_radius, 
                                     y₁ = bottom, y₂ = depth)
    elseif radius < 11
        return depth
    elseif radius < 12
        return cubic_profile(radius, x₁ = right_inner_radius, x₂ = right_outer_radius, 
                                     y₁ = depth, y₂ = bottom)
    else
        return bottom
    end
end

# Full bathymetry!
function bathymetry_with_ridge(x, y; params::NeverWorldBathymetryParameters; 
                               longitudinal_extent = 60, latitude = (-70, 70)) 
    if x < 5 || x > 55
        if x < 0 
           x = 0.0
        end
        if x > longitudinal_extent
           x = longitudinal_extent
        end

        channel_south = params.channel_south_edge
        channel_north = params.channel_north_edge
        bottom = params.bottom

        if y > channel_south && y < channel_north
            return  max(scotia_arc(x, y, params.scotia_arc, bottom), 
                        coastal_shelf_x(sqrt(x^2 + (y - channel_south)^2), params.shelves, bottom),
                        coastal_shelf_x(sqrt(x^2 + (y - channel_north)^2), params.shelves, bottom), 
                        coastal_shelf_x(sqrt((longitudinal_extent - x)^2 + (y - channel_south)^2), params.shelves, bottom),
                        coastal_shelf_x(sqrt((longitudinal_extent - x)^2 + (y - channel_north)^2), params.shelves, bottom))
        else
            return max(coastal_shelf_x(x, params.shelves, bottom), 
                       coastal_shelf_x(longitudinal_extent - x, params.shelves, bottom), 
                       coastal_shelf_y(-latitude[1] + y, params.shelves, bottom),
                       coastal_shelf_y(latitude[2]  - y, params.shelves, bottom),
                       bottom_ridge_xy(x, y, params.ridge, bottom), 
                       bottom_ridge_xy(longitudinal_extent - x, y, params.ridge, bottom), 
                       scotia_arc(x, y, params.scotia_arc, bottom))
        end
    else
        return max(coastal_shelf_x(x, params.shelves, bottom), 
                   coastal_shelf_x(longitudinal_extent - x, params.shelves, bottom),  
                   coastal_shelf_y(-latitude[1] + y, params.shelves, bottom), 
                   coastal_shelf_y(latitude[2]  - y, params.shelves, bottom), 
                   bottom_ridge_xy(x, y, params.ridge, bottom), 
                   bottom_ridge_xy(longitudinal_extent - x, y, params.ridge, bottom),  
                   scotia_arc(x, y, params.scotia_arc, bottom))
    end
end
