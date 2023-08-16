using Oceananigans
using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans.Fields: @compute

dir = "/storage2/WenoNeverworldData/"
prefix    = "weno_fourth_check"
variables = ("u", "v", "w", "b")
# stride = 3

# Load FieldTimeSeries of "u", "v", "w", "b" and compute integrated quantities
@info "loading files and calculating integrated variables"
neverworld_fields = all_fieldtimeseries(prefix, dir; checkpointer = true, variables, number_files = 20);
# kinetic_energy = Diagnostics.integral_kinetic_energy(neverworld_fields[:u], neverworld_fields[:b]; stride)
# integral_heat_content = Diagnostics.heat_content(neverworld_fields[:b]; stride)

grid  = neverworld_fields[:u].grid
times = neverworld_fields[:u].times
Nt = length(times)
iterations    = Nt-20:Nt
average_times = times[iterations]

u = neverworld_fields[:u]
v = neverworld_fields[:v]
w = neverworld_fields[:w]
b = neverworld_fields[:b]

# Average relevant variables
@info "Averaging relevant variables over" average_times
ū = Diagnostics.time_average(u, iterations)
v̄ = Diagnostics.time_average(v, iterations)
w_avg = Diagnostics.time_average(w, iterations)
b̄ = Diagnostics.time_average(b, iterations)
ψ = Diagnostics.calculate_eulerian_MOC(v̄)

# @compute b̄ᵢ = Field(Average(b̄, dims = 1))

loc = location(u^2 + v^2)
EKE = FieldTimeSeries{loc...}(grid, average_times)
IKE = FieldTimeSeries{loc...}(grid, average_times)

@compute MKE = Field(0.5 * (ū^2 + v̄^2))

u′ = FieldTimeSeries{location(u)...}(grid, average_times)
v′ = FieldTimeSeries{location(v)...}(grid, average_times)
w′ = FieldTimeSeries{location(w)...}(grid, average_times)
b′ = FieldTimeSeries{location(b)...}(grid, average_times)

@info "Computing fluctuating fields" 
for (idx, iter) in enumerate(iterations)
    @info "at iter $iter of $(iterations[end])"
    set!(u′[idx], u[iter] - ū)
    set!(v′[idx], v[iter] - v̄)
    set!(w′[idx], w[iter] - w_avg)
    set!(b′[idx], b[iter] - b̄)

    set!(EKE[idx], 0.5 * (u′[idx]^2 + v′[idx]^2))
    set!(IKE[idx], 0.5 * (u[idx]^2 + v[idx]^2))
end

@compute MKE_avg = Field(Average(MKE, dims = (1, 3)))
@compute EKE_avg = Field(Average(EKE, dims = (1, 3)))

v′b′ = FieldTimeSeries{location(v[1] * b[1])...}(grid, average_times)
w′b′ = FieldTimeSeries{location(w[1] * b[1])...}(grid, average_times)

@info "Computing eddy fluxes" 
for (idx, iter) in enumerate(iterations)
    @info "at iter $iter of $(iterations[end])"
    set!(v′b′[idx], v′[iter] * b′[iter])
    set!(w′b′[idx], v′[iter] * b′[iter])
end

@info "Computation of residual circulation"
vb_avg = Diagnostics.time_average(v′b′)
wb_avg = Diagnostics.time_average(w′b′)

@compute vbᵢ_avg = Field(Average(vb_avg, dims = 1))
@compute wbᵢ_avg = Field(Average(wb_avg, dims = 1))

@compute ∂ybᵢ = Field(∂y(b̄ᵢ)) 
@compute ∂zbᵢ = Field(∂z(b̄ᵢ)) 

@compute μ = Field(- wbᵢ_avg / vbᵢ_avg * ∂zbᵢ / ∂ybᵢ)