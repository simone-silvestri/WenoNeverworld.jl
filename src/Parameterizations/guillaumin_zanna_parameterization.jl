using Flux
using Flux: Conv, relu, Chain
using Flux.Optimise: softplus
using JLD2 
using OffsetArrays
using Adapt

using Oceananigans: architecture
import Oceananigans: on_architecture

@inline flux_architecture(::GPU) = Flux.gpu
@inline flux_architecture(::CPU) = Flux.cpu

struct NNbackscatteringClosure{NN, FT} <: AbstractTurbulenceClosure{ExplicitTimeDiscretization, 2}
    nn  :: NN # the convolutional neural network that computes `nn(u, v) -> (Su, Sv)`
    u★  :: FT # scaling constant for the zonal velocity
    v★  :: FT # scaling constant for the meridional velocity
    Su★ :: FT # scaling constant for the zonal subgrid scale forcing
    Sv★ :: FT # scaling constant for the meridional subgrid scale forcing
    min_value :: FT # Value to add to the activated elements (default: 0.0015).
    sampling :: Int
end

# On the GPU we throw away the NN since it is not intended to 
# be adapted for use inside kernels (and we do not need it inside kernels)
Adapt.adapt_structure(to, clo::NNbackscatteringClosure) = 
    NNbackscatteringClosure(nothing, clo.u★, clo.v★, clo.Su★, clo.Sv★, clo.min_value, clo.sampling)

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
                                 min_value = 0.0015,
                                 sampling = true)

    nn = getmodel(weight_path; architecture)
    
    u_scale  = convert(FT, u_scale)
    v_scale  = convert(FT, v_scale)
    Su_scale = convert(FT, Su_scale)
    Sv_scale = convert(FT, Sv_scale)
    min_value = convert(FT, min_value)

    return NNbackscatteringClosure(nn, u_scale,  v_scale, 
                                   Su_scale, Sv_scale, min_value, Int(sampling))
end

function DiffusivityFields(grid, tracer_names, bcs, ::NNbackscatteringClosure)
    arch = architecture(grid)

    # Inputs to the NN    
    utmp = CenterField(grid)

    # Outputs of the NN
    Su  = XFaceField(grid)
    Sv  = YFaceField(grid)

    Nx, Ny, Nz = size(utmp.data.parent)
    ox, oy, oz = utmp.data.offsets

    # Inpur work array -- 2 channels, where x, y, and z dimensions
    # are offset like the u and v fields, while the channel \
    # dimension is indexed from 1
    wrk_in = OffsetArray(zeros(Nx, Ny, 2, Nz), ox, oy, 0, oz)
    wrk_in = on_architecture(arch, wrk_in)

    # Output work array -- 4 channels, where x, y, and z dimensions
    # are offset like the u and v fields, while the channel \
    # dimension is indexed from 1
    wrk_out = OffsetArray(zeros(Nx, Ny, 4, Nz), ox, oy, 0, oz)
    wrk_out = on_architecture(arch, wrk_out)

    return (; Su, Sv, wrk_in, wrk_out)
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
    Su     = K.Su
    Sv     = K.Sv
    output = K.wrk_out

    # NN inputs
    input = K.wrk_in

    # Scaling parameters
    u★  = closure.u★
    v★  = closure.v★
    Su★ = closure.Su★
    Sv★ = closure.Sv★
    sampling = closure.sampling

    grid = u.grid
    arch = architecture(grid)

    launch!(arch, grid, parameters, _populate_input!, input, grid, u, v, u★, v★)

    #(w, h, 2, k) - Here we consider depth layers as batch as they are processed independently
    # Here we are allocating!!! (better to do inplace substitution if possible)
    # This step needs to be GPU-compatible, it's the last step we need to figure out
    output.parent .= closure.nn(input.parent) 

    # Apply the activation (the softplus function) pointwise
    launch!(arch, grid, parameters, _activation!, output, closure.min_value)

    # Sample the outputs on the correct device
    launch!(arch, grid, parameters, _sample_output!, Su, Sv, output, grid, Su★, Sv★, sampling)

    return nothing
end

# Interpolate velocities from staggered locations to centered locations
@kernel function _populate_input!(input, grid, u, v, u★, v★)
    i, j, k = @index(Global, NTuple)
    @inbounds input[i, j, 1, k] = ℑxᶜᵃᵃ(i, j, k, grid, u) * u★
    @inbounds input[i, j, 2, k] = ℑyᵃᶜᵃ(i, j, k, grid, v) * v★
end

@inline function sample_output_uᶜᶜᶜ(i, j, k, grid, output, Su★, sampling)
    @inbounds Su  = output[i, j, 1, k]
    @inbounds Spu = output[i, j, 3, k]
    return Su★ * (Su + sqrt(1 / Spu) * randn() * sampling)
end

@inline function sample_output_vᶜᶜᶜ(i, j, k, grid, output, Sv★, sampling)
    @inbounds Sv  = output[i, j, 2, k]
    @inbounds Spv = output[i, j, 4, k]
    return Sv★ * (Sv + sqrt(1 / Spv) * randn() * sampling)
end

# Compute the sampling output on centers and interpolate them onto the 
# staggered C-grid
@kernel function _sample_output!(Su, Sv, output, grid, Su★, Sv★, sampling)
    i, j, k = @index(Global, NTuple)

    @inbounds Su[i, j, k] = ℑxᶠᵃᵃ(i, j, k, grid, sample_output_uᶜᶜᶜ, output, Su★, sampling)
    @inbounds Sv[i, j, k] = ℑyᵃᶠᵃ(i, j, k, grid, sample_output_vᶜᶜᶜ, output, Sv★, sampling)
end

# Forcing in the u- and v- equations
@inline ∂ⱼ_τ₁ⱼ(i, j, k, grid, ::NNbackscatteringClosure, K, args...) = @inbounds K.Su[i, j, k]
@inline ∂ⱼ_τ₂ⱼ(i, j, k, grid, ::NNbackscatteringClosure, K, args...) = @inbounds K.Sv[i, j, k]
    
# No forcing term in the w-equation or in the tracer equations!
@inline ∂ⱼ_τ₃ⱼ(i, j, k, grid, ::NNbackscatteringClosure, args...)   = zero(grid)
@inline ∇_dot_qᶜ(i, j, k, grid, ::NNbackscatteringClosure, args...) = zero(grid)

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
@kernel function _activation!(x, min_value) 
    i, j, k = @index(Global, NTuple)
    @inbounds x[i, j, 3, k] = softplus(x[i, j, 3, k]) + min_value
    @inbounds x[i, j, 4, k] = softplus(x[i, j, 4, k]) + min_value
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
function getmodel(weight_path=nothing; architecture = CPU(), padding="init_zeros") 
    # Define the network structure

    if padding == "same"
        p5 = 2
        p3 = 1
        p_init = p5
    elseif padding == "init_zeros" 
        p5 = 0
        p3 = 0 
        p_init = 10
    else
        error("padding option $padding unknown")
    end


    model = Chain(
        Conv((5, 5), 2   => 128, pad = (p_init, p_init)), relu,
        Conv((5, 5), 128 => 64,  pad = (p5, p5)), relu,
        Conv((3, 3), 64  => 32,  pad = (p3, p3)), relu,
        Conv((3, 3), 32  => 32,  pad = (p3, p3)), relu,
        Conv((3, 3), 32  => 32,  pad = (p3, p3)), relu,
        Conv((3, 3), 32  => 32,  pad = (p3, p3)), relu,
        Conv((3, 3), 32  => 32,  pad = (p3, p3)), relu,
        Conv((3, 3), 32  => 4,   pad = (p3, p3)))
    if !isnothing(weight_path) 
        @info "Loading model : ", weight_path
        model_state = JLD2.load(weight_path, "model_state");
        Flux.loadmodel!(model, model_state);
    end
    flux_arch = flux_architecture(architecture)

    return model |> flux_arch
end
