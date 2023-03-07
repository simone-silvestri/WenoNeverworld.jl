# WenoNeverworld.jl

Modified neverworld2 setup from [(original neverworld2 setup)](https://egusphere.copernicus.org/preprints/2022/egusphere-2022-186/egusphere-2022-186.pdf) to 
1. study mesh convergence  
2. investigate numerics
3. calibrate and assess mesoscale parameterizations

## Example

Build and run a neverworld simulation 

```
using Oceananigans
using Oceananigans.Units
using WenoNeverworld

output_dir    = joinpath(@__DIR__, "../output_directory")
@show output_prefix = output_dir * "/neverworld_simulation_quarter_degree"

arch       = GPU()
resolution = 1/4

# Construct a grid without an atlantic ridge with λ ∈ (-2, 62), φ ∈ (-70, 0) and 69 exponential levels
grid = NeverworldGrid(arch, resolution)

# Start simulation from scratch
interp_init = false
init_file   = nothing

# Simulation parameters
Δt        = 10minutes
stop_time = 7000days

# Construct the neverworld simulation
simulation = weno_neverworld_simulation(; grid, Δt, stop_time, interp_init, init_file)

# Increase simulation Δt to 16 minutes after 300 days
increase_simulation_Δt!(simulation, cutoff_time = 100days,  new_Δt = 12minutes)
increase_simulation_Δt!(simulation, cutoff_time = 200days,  new_Δt = 14minutes)
increase_simulation_Δt!(simulation, cutoff_time = 300days,  new_Δt = 16minutes)

# Attach standard outputs to the simulation
standard_outputs!(simulation, output_prefix)

# Let's run!
run_simulation!(simulation; interp_init, init_file)
```

## Speed at the surface (1/8th WENO, 1/4er WENO, 1/4er Leith)

<img width="1067" alt="image" src="https://user-images.githubusercontent.com/33547697/223501084-cdc3b3e2-4ef2-4b4a-8db2-8f93fdb1133c.png">



