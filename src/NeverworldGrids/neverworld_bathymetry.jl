# The bathymetry is defined for a latitude range of -70 ≤ φ ≤ 70
# and a longitude range of 0 ≤ λ ≤ 60

# All quantities in the horizontal direction are specified in degrees and in the vertical in meters
Base.@kwdef struct ShelfParameters
    coast_length::Float64 = 0.5
    side_length::Float64  = 2.5
    length::Float64       = 2.5
    depth::Float64        = 200
end

Base.@kwdef struct RidgeParameters
    side_length::Float64    = 9
    top_length::Float64     = 2
    longitude::Float64      = 30
    south_latitude::Float64 = -30
    slope_length::Float64   = 20
    depth::Float64          = 2000
end

Base.@kwdef struct ScotiaArcParameters
    left_inner_radius::Float64  = 8
    left_outer_radius::Float64  = 9
    right_inner_radius::Float64 = 11
    right_outer_radius::Float64 = 12
    center_latitude::Float64    = 50
    depth::Float64              = 2000
end

Base.@kwdef struct NeverWorldBathymetryParameters
    shelves                     = ShelfParameters()
    scotia_arc                  = ScotiaArcParameters()
    channel_south_edge::Float64 = - 59
    channel_north_edge::Float64 = - 41
    bottom::Float64             = - 4000
    ridge                       = nothing
end
    
## define the coasts
function coastal_shelf_x(x, params, bottom) 

    coast       = params.coast_length
    length      = params.length
    side_length = params.side_length + length
    depth       = - params.depth

    if x < coast
        return 0.0
    elseif x < length
        return depth
    elseif x < side_length
        return cubic_interpolate(x, x₁ = length, x₂ = side_length, y₁ = depth, y₂ = bottom)
    else        
        return bottom
    end
end

function sharp_coast_x(x, params, bottom) 

    coast  = params.coast_length

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
    depth       = - params.depth

    if x < coast
        return cubic_interpolate(x, x₁ = 0.0, x₂ = coast, y₁ = 0.0, y₂ = depth)
    elseif x < length
        return depth
    elseif x < side_length
        return cubic_interpolate(x, x₁ = length, x₂ = side_length, y₁ = depth, y₂ = bottom)
    else        
        return bottom
    end
end

# Bottom ridge
function bottom_ridge_x(x, params, bottom)

    center     = params.longitude
    
    top_left   = center - params.top_length/2
    top_right  = center + params.top_length/2

    bot_left  = center - params.top_length/2 - params.side_length
    bot_right = center + params.top_length/2 + params.side_length

    depth = - params.depth

    if x < bot_left
        return bottom
    elseif x < top_left
        return cubic_interpolate(x, x₁ = bot_left, x₂ = top_left, y₁ = bottom, y₂ = depth)
    elseif x < top_right
        return depth
    elseif x < bot_right
        return cubic_interpolate(x, x₁ = top_right, x₂ = bot_right, y₁ = depth, y₂ = bottom)
    else 
        return bottom
    end
end
       
""" smoothed coasts for the inlet and outlet of the channel """
bottom_ridge_xy(x, y, ::Nothing, bottom) = bottom

function bottom_ridge_xy(x, y, params, bottom)
    sl = params.south_latitude
    bl = params.south_latitude - params.slope_length

    if y > sl
        return bottom_ridge_x(x, params, bottom)
    elseif y > bl
        return cubic_interpolate(y, x₁ = sl, x₂ = bl, y₁ = bottom_ridge_x(x, params, bottom), y₂ = bottom)
    else
        return bottom
    end
end

scotia_arc(x, y, ::Nothing, bottom) = bottom

# Scotia arc
function scotia_arc(x, y, params, bottom)
    
    left_inner_radius  = params.left_inner_radius
    left_outer_radius  = params.left_outer_radius
    right_inner_radius = params.right_inner_radius
    right_outer_radius = params.right_outer_radius
    mid_point = params.center_latitude
    depth = - params.depth

    radius = sqrt(x^2 + (y + mid_point)^2)
    if radius < left_inner_radius
        return bottom
    elseif radius < left_outer_radius
        return cubic_interpolate(radius, x₁ = left_inner_radius, x₂ = left_outer_radius, 
                                         y₁ = bottom, y₂ = depth)
    elseif radius < right_inner_radius
        return depth
    elseif radius < right_outer_radius
        return cubic_interpolate(radius, x₁ = right_inner_radius, x₂ = right_outer_radius, 
                                         y₁ = depth, y₂ = bottom)
    else
        return bottom
    end
end

# Full bathymetry!
function neverworld_bathymetry(x, y, params::NeverWorldBathymetryParameters; 
                               longitudinal_extent = 60, latitude = (-70, 70)) 
    
    channel_south = params.channel_south_edge
    channel_north = params.channel_north_edge
    bottom = params.bottom
                           
    if x < 5 || x > 55
        if x < 0 
           x = 0.0
        end
        if x > longitudinal_extent
           x = longitudinal_extent
        end

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
