using KernelAbstractions: @kernel, @index
using Oceananigans.Grids: architecture
using Oceananigans.Utils
using Oceananigans.BoundaryConditions

# Maybe we can remove this propagate field in lieu of a diffusion, 
# Still we'll need to do this a couple of steps on the original grid
@kernel function _propagate_field!(field, tmp_field)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        nw = field[i - 1, j, k]
        ns = field[i, j - 1, k]
        ne = field[i + 1, j, k]
        nn = field[i, j + 1, k]
        nb = (nw, ne, nn, ns)
    end

    counter = 0
    cumsum  = zero(eltype(field))

    for n in nb
        counter += ifelse(isnan(n), 0, 1)
        cumsum  += ifelse(isnan(n), 0, n)
    end

    @inbounds tmp_field[i, j, k] = ifelse(cumsum == 0, NaN, cumsum / counter)
end

@kernel function _substitute_values!(field, tmp_field)
    i, j, k = @index(Global, NTuple)
    @inbounds needs_inpainting = isnan(field[i, j, k])
    @inbounds field[i, j, k] = ifelse(needs_inpainting, tmp_field[i, j, k], field[i, j, k])
end

@kernel function _nan_field!(field)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] = ifelse(field[i, j, k] == 0, NaN, field[i, j, k])
end

propagate_horizontally!(field, ::Nothing, tmp_field=deepcopy(field); kw...) = field

function propagating(field, iter, maxiter)
    mask_sum = sum(field)
    return isnan(mask_sum) && iter < maxiter
end

@kernel function _fill_nans!(field)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] = ifelse(isnan(field[i, j, k]), 0, field[i, j, k])
end

""" 
    propagate_horizontally!(field, mask [, tmp_field=deepcopy(field)]; max_iter = Inf)

Horizontally propagate the values of `field` into the `mask`.
In other words, cells where `mask[i, j, k] == false` are preserved,
and cells where `mask[i, j, k] == true` are painted over.
"""
function propagate_horizontally!(field, tmp_field=deepcopy(field); maxiter = Inf) 
    iter  = 0
    grid  = field.grid
    arch  = architecture(grid)
    
    launch!(arch, grid, :xyz, _nan_field!, field)
    fill_halo_regions!(field)

    # Need temporary field to avoid a race condition
    parent(tmp_field) .= parent(field)

    while propagating(field, iter, maxiter)
        launch!(arch, grid, :xyz, _propagate_field!,   field, tmp_field)
        launch!(arch, grid, :xyz, _substitute_values!, field, tmp_field)
        iter += 1
        @debug "Propagate pass $iter with sum $(sum(parent(field)))"
    end

    launch!(arch, grid, :xyz, _fill_nans!, field)

    fill_halo_regions!(field)

    return field
end