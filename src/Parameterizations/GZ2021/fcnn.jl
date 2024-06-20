using Flux: Conv, relu, Chain
using Flux.Optimise: softplus

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
    get_model(weight_path=nothing)

Defines and optionally loads the weights for network from 

    Guillaumin, A. P., & Zanna, L. (2021). Stochastic deep learning parameterization of ocean momentum forcing. Journal of Advances in Modeling Earth Systems, 13(9), e2021MS002534.
    https://github.com/chzhangudel/Forpy_CNN_GZ21

# Arguments
- `weight_path`: Path to the file containing the model weights (default: nothing).

# Returns
- The constructed model, with weights loaded if `weight_path` is provided.
"""
function get_model(weight_path=nothing) 
    # Define the network structure
    model = Chain(
        Conv((5, 5), 2 => 128, pad = (2, 2)), relu,
        Conv((5, 5), 128 => 64, pad = (2, 2)), relu,
        Conv((3, 3), 64 => 32, pad = (1, 1)), relu,
        Conv((3, 3), 32 => 32, pad = (1, 1)), relu,
        Conv((3, 3), 32 => 32, pad = (1, 1)), relu,
        Conv((3, 3), 32 => 32, pad = (1, 1)), relu,
        Conv((3, 3), 32 => 32, pad = (1, 1)), relu,
        Conv((3, 3), 32 => 4, pad = (1, 1)))
    if !isnothing(weight_path) 
        @info "Loading model : ", weight_path
        model_state = JLD2.load(weight_path, "model_state");
        Flux.loadmodel!(model, model_state);
    end
    return model
end


"""
    subgrid_forcing(u, v; u_scale=10.0, v_scale=10.0, Su_scale=1e-7, Sv_scale=1e-7, sampling=true)

Computes the subgrid forcing for zonal and meridional velocities using a neural network model.
Values from scale : https://github.com/chzhangudel/Forpy_CNN_GZ21/blob/smartsim/testNN.py

```julia
model=get_model("model_weights.jld2")
input = randn(Float32, 200, 128, 2, 10); # i, j, c (u; v), k (depth levels)
out = activation(model(input));
````

# Arguments
- `u`: Zonal velocity tensor (i, j, k (depth levels)).
- `v`: Meridional velocity tensor (i, j, k (depth levels)).
- `u_scale`: Scaling factor for the zonal velocity (default: 10.0).
- `v_scale`: Scaling factor for the meridional velocity (default: 10.0).
- `Su_scale`: Scaling factor for the zonal subgrid forcing (default: 1e-7).
- `Sv_scale`: Scaling factor for the meridional subgrid forcing (default: 1e-7).
- `sampling`: Boolean indicating if noise should be added to the subgrid forcing (default: true).


# Returns
- `Su`: Zonal subgrid forcing tensor (i, j, k (depth levels)).
- `Sv`: Meridional subgrid forcing tensor (i, j, k (depth levels)).

"""
function subgrid_forcing(u, v; u_scale=10.0, v_scale=10.0, Su_scale=1e-7, Sv_scale=1e-7, sampling=true)
    out = model(stack([u .* u_scale, v .* v_scale], dims=3)) #(w, h, 2, k) - Here we consider depth layers as batch as they are processed independently
    out = activation(out)
    Su, Sv, Spu, Spv = out[:,:,1,:], out[:,:,2,:], out[:,:,3,:], out[:,:,4,:]
    Su = Su_scale * ( Su + sqrt.(1 ./Spu).*randn(size(Spu)) * sampling )
    Sv = Sv_scale * ( Sv + sqrt.(1 ./Spv).*randn(size(Spv)) * sampling )
    return Su, Sv
end



