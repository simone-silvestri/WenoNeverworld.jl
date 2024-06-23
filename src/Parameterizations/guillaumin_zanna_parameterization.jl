using Flux
using Flux: Conv, relu, Chain
using Flux.Optimise: softplus
using JLD2 

using Oceananigans: architecture
import Oceananigans: on_architecture

import Oceananigans.TurbulenceClosures: 
                        ∂ⱼ_τ₁ⱼ, 
                        ∂ⱼ_τ₂ⱼ, 
                        ∂ⱼ_τ₃ⱼ,
                        ∇_dot_qᶜ

import Oceananigans.TurbulenceClosures: compute_diffusivities!, DiffusivityFields

struct NNbackscatteringClosure{NN, FT} <: AbstractTurbulenceClosure{ExplicitTimeDiscretization, 2}
    nn  :: NN # the convolutional neural network that computes `nn(u, v) -> (Su, Sv)`
    u★  :: FT # scaling constant for the zonal velocity
    v★  :: FT # scaling constant for the meridional velocity
    Su★ :: FT # scaling constant for the zonal subgrid scale forcing
    Sv★ :: FT # scaling constant for the meridional subgrid scale forcing
    sampling :: Int
end

"""
    NNbackscatteringClosure(FT::DataType = Float64; 
                        weight_path = nothing,
                        u_scale = 10,
                        v_scale = 10,
                        Su_scale = 1e-7, 
                        Sv_scale = 1e-7, 
                        sampling = true)

Constructs a subgrid-scale closure implemented as a neural network. This closure computes the 
subgrid-scale forcing in the `compute_diffusivities!` step and then applies them to the `u-momentum` 
and `v-momentum` equations extending the flux divergence functions `∂ⱼ_τ₁ⱼ` and `∂ⱼ_τ₂ⱼ`

# Arguments
============
- `FT::DataType`: The data type to use for the model parameters. Defaults to `Float64`.

# Keyword Arguments
===================
- `weight_path`: The path to the pre-trained weights of the neural network model. Defaults to `nothing`.
- `u_scale`: The scaling factor for the u-component of the velocity field. Defaults to `10`.
- `v_scale`: The scaling factor for the v-component of the velocity field. Defaults to `10`.
- `Su_scale`: The scaling factor for the subgrid-scale u-component forcing. Defaults to `1e-7`.
- `Sv_scale`: The scaling factor for the subgrid-scale v-component forcing. Defaults to `1e-7`.
- `sampling`: A boolean indicating whether to use sampling during the forward pass of the neural network. Defaults to `true`.
"""
function NNbackscatteringClosure(FT::DataType = Float64; 
                                 architecture = CPU(),
                                 weight_path = nothing,
                                 u_scale = 10,
                                 v_scale = 10,
                                 Su_scale = 1e-7, 
                                 Sv_scale = 1e-7, 
                                 sampling = true)

    nn = getmodel(weight_path; architecture)
    
    u_scale  = convert(FT, u_scale)
    v_scale  = convert(FT, v_scale)
    Su_scale = convert(FT, Su_scale)
    Sv_scale = convert(FT, Sv_scale)

    return NNbackscatteringClosure(nn, u_scale,  v_scale, 
                                   Su_scale, Sv_scale, Int(sampling))
end

function DiffusivityFields(grid, tracer_names, bcs, ::NNbackscatteringClosure)
    arch = architecture(grid)

    # Inputs to the NN    
    uᶜᶜᶜ = CenterField(grid)
    vᶜᶜᶜ = CenterField(grid)

    # Outputs of the NN
    Su  = XFaceField(grid)
    Sv  = YFaceField(grid)

    # # Work array (4 channels)
    # Nx, Ny, Nz = size(grid)
    # wrk = zeros(Nx, Ny, 4, Nz)
    # wrk = on_architecture(arch, wrk)

    return (; uᶜᶜᶜ, vᶜᶜᶜ, Su, Sv) #, wrk)
end

#####
##### Forcing-specific functions 
#####

"""
    compute_diffusivities!(K, closure::NNbackscatteringClosure, model; parameters = :xyz)

Computes the subgrid forcing for zonal and meridional velocities using a neural network model.
Values from scale : https://github.com/chzhangudel/Forpy_CNN_GZ21/blob/smartsim/testNN.py
"""
# Calculate forcing terms during the `compute_diffusivities!` step (before calculating tendencies)
function compute_diffusivities!(K, closure::NNbackscatteringClosure, model; parameters = :xyz)
    arch = model.architecture
    grid = model.grid
    u, v, _ = model.velocities

    # NN outputs
    Su = K.Su
    Sv = K.Sv

    # NN inputs
    uᶜᶜᶜ = K.uᶜᶜᶜ
    vᶜᶜᶜ = K.vᶜᶜᶜ

    # Scaling parameters
    u★  = closure.u★
    v★  = closure.v★
    Su★ = closure.Su★
    Sv★ = closure.Sv★
    sampling = closure.sampling

    grid = u.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz, _scaled_center_velocities!, uᶜᶜᶜ, vᶜᶜᶜ, grid, u, v, u★, v★)

    #(w, h, 2, k) - Here we consider depth layers as batch as they are processed independently
    # Here we are allocating!!! (better to do inplace substitution if possible)
    out = closure.nn(stack([interior(uᶜᶜᶜ), interior(vᶜᶜᶜ)], dims=3)) 
    out = activation(out)
    
    # Sample the outputs on the correct device
    launch!(arch, grid, parameters, _sample_output!, Su, Sv, out, grid, Su★, Sv★, sampling)

    return nothing
end

# Interpolate velocities from staggered locations to centered locations
@kernel function _scaled_center_velocities!(uᶜᶜᶜ, vᶜᶜᶜ, grid, u, v, u★, v★)
    i, j, k = @index(Global, NTuple)
    @inbounds uᶜᶜᶜ[i, j, k] = ℑxᶜᵃᵃ(i, j, k, grid, u) * u★
    @inbounds vᶜᶜᶜ[i, j, k] = ℑyᵃᶜᵃ(i, j, k, grid, v) * v★
end

@inline function sample_output_uᶜᶜᶜ(i, j, k, grid, out, Su★, sampling)
    @inbounds Su  = out[i, j, 1, k]
    @inbounds Spu = out[i, j, 3, k]
    return Su★ * (Su + sqrt(1 / Spu) * randn() * sampling)
end

@inline function sample_output_vᶜᶜᶜ(i, j, k, grid, out, Sv★, sampling)
    @inbounds Sv  = out[i, j, 2, k]
    @inbounds Spv = out[i, j, 4, k]
    return Sv★ * (Sv + sqrt(1 / Spv) * randn() * sampling)
end

# Compute the sampling output on centers and interpolate them onto the 
# staggered C-grid
@kernel function _sample_output!(Su, Sv, out, grid, Su★, Sv★, sampling)
    i, j, k = @index(Global, NTuple)

    @inbounds Su[i, j, k] = ℑxᶠᵃᵃ(i, j, k, grid, sample_output_uᶜᶜᶜ, out, Su★, sampling)
    @inbounds Sv[i, j, k] = ℑyᵃᶠᵃ(i, j, k, grid, sample_output_vᶜᶜᶜ, out, Sv★, sampling)
end

# Forcing in the u- and v- equations
@inline ∂ⱼ_τ₁ⱼ(i, j, k, grid, closure::NNbackscatteringClosure, K, args...) = @inbounds K.Su[i, j, k]
@inline ∂ⱼ_τ₂ⱼ(i, j, k, grid, closure::NNbackscatteringClosure, K, args...) = @inbounds K.Sv[i, j, k]
    
# No forcing term in the w-equation or in the tracer equations!
@inline ∂ⱼ_τ₃ⱼ(i, j, k, grid, closure::NNbackscatteringClosure, args...)   = zero(grid)
@inline ∇_dot_qᶜ(i, j, k, grid, closure::NNbackscatteringClosure, args...) = zero(grid)

#####
##### NN-specific functions
#####

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
    view(out, :, :, precision_indices, :) .= softplus.(view(x, :, :, precision_indices, :)) .+ min_value
    return out
end

"""
    getmodel(weight_path=nothing)

Defines and optionally loads the weights for network from 

    Guillaumin, A. P., & Zanna, L. (2021). Stochastic deep learning parameterization of ocean momentum forcing. Journal of Advances in Modeling Earth Systems, 13(9), e2021MS002534.
    https://github.com/chzhangudel/Forpy_CNN_GZ21

# Arguments
============
- `weight_path`: Path to the file containing the model weights (default: nothing).

# Keyword Arguments
===================
- `architecture`: the architecture on which the model runs: either `CPU()` or `GPU()`

# Returns
- The constructed model, with weights loaded if `weight_path` is provided.
"""
function getmodel(weight_path=nothing; architecture = CPU()) 
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
    return on_architecture(architecture, model)
end

# Make the `Chain` structure GPU-compatible by converting all the
# concrete arrays and data structures to their GPU-compatible counterparts
# In this case, we only need to convert the `weight`s and the `bias`es
# TODO: make sure there is no better way to do this step already implemented in `Flux`
function on_architecture(arch, nn :: Chain)
    new_layers = []

    for layer in nn.layers
        if layer isa Function 
            push!(new_layers, layer)
        else
            weight = on_architecture(arch, layer.weight)
            bias   = on_architecture(arch, layer.bias)

            new_layer = Conv(weight, bias, layer.σ; stride = layer.stride, 
                                                       pad = layer.pad, 
                                                  dilation = layer.dilation)

            push!(new_layers, new_layer)
        end
    end

    return Chain(tuple(new_layers...))
end