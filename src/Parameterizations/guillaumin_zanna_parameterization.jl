using Flux
using Flux: Conv, relu, Chain
using Flux.Optimise: softplus
using JLD2 

import Oceananigans.TurbulenceClosures: 
                        ∂ⱼ_τ₁ⱼ, 
                        ∂ⱼ_τ₂ⱼ, 
                        ∂ⱼ_τ₃ⱼ,
                        ∇_dot_qᶜ

import Oceananigans.TurbulenceClosures: compute_diffusivities!, DiffusivityFields

struct NNSubgridSaleForcing{NN, FT} <: AbstractTurbulenceClosure{ExplicitTimeDiscretization, 2}
    nn :: NN
    u_scale  :: FT
    v_scale  :: FT
    Su_scale :: FT
    Sv_scale :: FT
    sampling :: Int
end

function NNSubgridSaleForcing(FT::DataType = Float64; 
                              weight_path = nothing,
                              u_scale = 10,
                              v_scale = 10,
                              Su_scale = 1e-7, 
                              Sv_scale = 1e-7, 
                              sampling = true)

    nn = getmodel(weight_path)
    
    u_scale  = convert(FT, u_scale)
    v_scale  = convert(FT, v_scale)
    Su_scale = convert(FT, Su_scale)
    Sv_scale = convert(FT, Sv_scale)

    return NNSubgridSaleForcing(nn, 
                                u_scale,  v_scale, 
                                Su_scale, Sv_scale, Int(sampling))
end

DiffusivityFields(grid, tracer_names, bcs, ::NNSubgridSaleForcing) = 
                (; Su  = XFaceField(grid),
                   Sv  = YFaceField(grid))


"""
    compute_diffusivities!(K, closure::NNSubgridSaleForcing, model; parameters = :xyz)

Computes the subgrid forcing for zonal and meridional velocities using a neural network model.
Values from scale : https://github.com/chzhangudel/Forpy_CNN_GZ21/blob/smartsim/testNN.py
"""
# Calculate forcing terms during the `compute_diffusivities!` step (before calculating tendencies)
function compute_diffusivities!(K, closure::NNSubgridSaleForcing, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    u, v, _ = model.velocities

    grid = u.grid
    arch = architecture(grid)

    # Forcing fields
    Su = K.Su
    Sv = K.Sv

    # Scaling parameters
    u_scale  = closure.u_scale
    v_scale  = closure.v_scale
    Su_scale = closure.Su_scale
    Sv_scale = closure.Sv_scale

    #(w, h, 2, k) - Here we consider depth layers as batch as they are processed independently
    # Here we are allocating!!! (better to do inplace substitution if possible)
    out = closure.nn(stack([u .* u_scale, v .* v_scale], dims=3)) 
    out = activation(out)
    
    # Sample the outputs on the correct device
    launch!(arch, grid, parameters, _sample_output!, Su, Sv, out, Su_scale, Sv_scale, sampling)

    return nothing
end

@kernel function _sample_output!(Su, Sv, out, Su_scale, Sv_scale, sampling)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        Su[i, j, k] = out[i, j, 1, k]
        Sv[i, j, k] = out[i, j, 2, k]
        
        Spu = out[i, j, 3, k]
        Spv = out[i, j, 4, k]

        Su[i, j, k] = Su_scale * (Su[i, j, k] + sqrt(1 / Spu) * randn(Spu) * sampling)
        Sv[i, j, k] = Sv_scale * (Sv[i, j, k] + sqrt(1 / Spv) * randn(Spv) * sampling)
    end
end

# Forcing in the u- and v- equations
@inline ∂ⱼ_τ₁ⱼ(i, j, k, grid, closure::NNSubgridSaleForcing, K, args...) = @inbounds K.Su[i, j, k]
@inline ∂ⱼ_τ₂ⱼ(i, j, k, grid, closure::NNSubgridSaleForcing, K, args...) = @inbounds K.Sv[i, j, k]
    
# No forcing term in the w-equation or in the tracer equations!
@inline ∂ⱼ_τ₃ⱼ(i, j, k, grid, closure::NNSubgridSaleForcing, args...)   = zero(grid)
@inline ∇_dot_qᶜ(i, j, k, grid, closure::NNSubgridSaleForcing, args...) = zero(grid)


"""
    activation(x; precision_indices=3:4, min_value=0.0015)

Applies the softplus activation function to the specified indices of the input tensor `x`, and adds a minimum value to them.

# Arguments
- `x`: Input tensor. (i, j, c (channels), k (depth levels)) 
- `precision_indices`: Indices along the third dimension to which the activation function is applied (default: 3:4).
- `min_value`: Value to add to the activated elements (default: 0.0015).

# Returns
- A new tensor with the activation applied to the specified indices.
"""
function activation(x; precision_indices=3:4, min_value=0.0015)
    out = copy(x) # If we want to avoid inplace modification
    out[:, :, precision_indices, :] .= softplus.(x[:, :, precision_indices, :]) .+ min_value
    return out
end

"""
    getmodel(weight_path=nothing)

Defines and optionally loads the weights for network from 

    Guillaumin, A. P., & Zanna, L. (2021). Stochastic deep learning parameterization of ocean momentum forcing. Journal of Advances in Modeling Earth Systems, 13(9), e2021MS002534.
    https://github.com/chzhangudel/Forpy_CNN_GZ21

# Arguments
- `weight_path`: Path to the file containing the model weights (default: nothing).

# Returns
- The constructed model, with weights loaded if `weight_path` is provided.
"""
function getmodel(weight_path=nothing) 
    # Define the network structure
    model = Chain(
        Conv((5, 5), 2   => 128, pad = (2, 2)), relu,
        Conv((5, 5), 128 => 64,  pad = (2, 2)), relu,
        Conv((3, 3), 64  => 32,  pad = (1, 1)), relu,
        Conv((3, 3), 32  => 32,  pad = (1, 1)), relu,
        Conv((3, 3), 32  => 32,  pad = (1, 1)), relu,
        Conv((3, 3), 32  => 32,  pad = (1, 1)), relu,
        Conv((3, 3), 32  => 32,  pad = (1, 1)), relu,
        Conv((3, 3), 32  => 4,   pad = (1, 1)))
    if !isnothing(weight_path) 
        @info "Loading model : ", weight_path
        model_state = JLD2.load(weight_path, "model_state");
        Flux.loadmodel!(model, model_state);
    end
    return model
end