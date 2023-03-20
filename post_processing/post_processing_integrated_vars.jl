using WenoNeverworld
using WenoNeverworld.Diagnostics
using Oceananigans
using Oceananigans.Units
using Oceananigans.Operators
using GLMakie
using Statistics: mean

"""
Time-averaged variables
"""

stride    = 10
variables = ("u", "v", "u2", "v2", "b", "wb", "w")

weno_avg   = all_fieldtimeseries("UpwindVectorInvariantSchemes/quarter"; variables)
center_avg = all_fieldtimeseries("UpwindVectorInvariantSchemes/quarter_centered"; variables)
leith_avg  = all_fieldtimeseries("UpwindVectorInvariantSchemes/quarter_leith"; variables)

Nt = length(center_avg[:u].times) 

@info "computing integrated heat content..."
center_heat = WenoNeverworld.Diagnostics.heat_content(center_avg[:b]; stride)
weno_heat   = WenoNeverworld.Diagnostics.heat_content(weno_avg[:b];   stride)

center_acc = WenoNeverworld.Diagnostics.ACC_transport(center_avg[:u]; stride)
weno_acc   = WenoNeverworld.Diagnostics.ACC_transport(weno_avg[:u];   stride)

Nt = Int(length(center_avg[:u].times) ÷ stride)

using WenoNeverworld.Diagnostics: DensityField

@info "computing equilibrium height z★..."
grid = center_avg[:u].grid

zc = Field((Center, Center, Center), grid)
zw = Field((Center, Center, Center), grid)

εc = Field((Center, Center, Center), grid)
εw = Field((Center, Center, Center), grid)
αc = Field((Center, Center, Center), grid)
αw = Field((Center, Center, Center), grid)

zfield = WenoNeverworld.Diagnostics.HeightField(grid)
vol        = VolumeField(grid)
area       = AreaField(grid, (Center, Center, Center))
total_area = sum(AreaField(grid))

RPEc = zeros(Nt)
RPEw = zeros(Nt)

APEc = zeros(Nt)
APEw = zeros(Nt)

KEc = zeros(Nt)
KEw = zeros(Nt)

"""
Eₖ and Eₐ budgets (only the important ones)
"""

ϕzc = zeros(Nt)
ϕzw = zeros(Nt)

ϕb1c = zeros(Nt)
ϕb1w = zeros(Nt)

ϕb2c = zeros(Nt)
ϕb2w = zeros(Nt)

ϕic = zeros(Nt)
ϕiw = zeros(Nt)

ρ₀ = 1000.0
g  = 9.80655

@info "computing resting and available potential energy density..."
for t in 1:Nt
    @info "time $t of $(Nt)"

    u2c = center_avg[:u2][(t - 1) * stride + 1]
    v2c = center_avg[:v2][(t - 1) * stride + 1]
    
    u2w = weno_avg[:u2][(t - 1) * stride + 1]
    v2w = weno_avg[:v2][(t - 1) * stride + 1]

    wc  = center_avg[:w][(t - 1) * stride + 1]
    wbc = center_avg[:wb][(t - 1) * stride + 1]
    
    ww  = weno_avg[:w][(t - 1) * stride + 1]
    wbw = weno_avg[:wb][(t - 1) * stride + 1]

    bc = center_avg[:b][(t - 1) * stride + 1]
    bw = weno_avg[:b][(t - 1) * stride + 1]

    ρc = DensityField(bc; ρ₀, g)
    ρw = DensityField(bw; ρ₀, g)

    @info "computing integrated energetics..."
    WenoNeverworld.Diagnostics.calculate_z★!(zc, bc, vol, total_area)
    WenoNeverworld.Diagnostics.calculate_z★!(zw, bw, vol, total_area)

    set!(εc, compute!(Field(zc * ρc)))
    set!(εw, compute!(Field(zw * ρw)))

    set!(αc, compute!(Field((- zfield - zc) * ρc)))
    set!(αw, compute!(Field((- zfield - zw) * ρc)))

    RPEc[t] = g * sum(compute!(Field(εc * vol))) 
    RPEw[t] = g * sum(compute!(Field(εw * vol))) 

    APEc[t] = g * sum(compute!(Field(αc * vol)))
    APEw[t] = g * sum(compute!(Field(αw * vol)))

    KEc[t] = ρ₀ * sum(compute!(Field(0.5 * (u2c + v2c) * vol)))
    KEw[t] = ρ₀ * sum(compute!(Field(0.5 * (u2w + v2w) * vol)))

    @info "computing budgets..."
    ∂zρc = compute!(Field(∂z(ρc)))
    ∂zρw = compute!(Field(∂z(ρw)))

    field_to_integ1c = compute!(Field(- zfield * ∂zρc * vol))
    field_to_integ1w = compute!(Field(- zfield * ∂zρw * vol))

    field_to_integ2c = compute!(Field(zc * ∂zρc * vol))
    field_to_integ2w = compute!(Field(zw * ∂zρw * vol))

    integ1c = sum(field_to_integ1c, dims = (1, 2))
    integ1w = sum(field_to_integ1w, dims = (1, 2))
    
    integ2c = sum(field_to_integ2c, dims = (1, 2))
    integ2w = sum(field_to_integ2w, dims = (1, 2))

    ϕb1c[t] = 1e-5 * g * (integ1c[1, 1, 1] - integ1c[1, 1, 69])
    ϕb1w[t] = 1e-5 * g * (integ1w[1, 1, 1] - integ1w[1, 1, 69])
    
    ϕb2c[t] = - 1e-5 * g * (integ2c[1, 1, 1] - integ2c[1, 1, 69])
    ϕb2w[t] = - 1e-5 * g * (integ2w[1, 1, 1] - integ2w[1, 1, 69])

    ρmeanc = mean(ρc, dims = (1, 2))
    ρmeanw = mean(ρw, dims = (1, 2))

    ϕic[t] = - 1e-5 * g * (ρmeanc[1, 1, 1] - ρmeanc[1, 1, 69])
    ϕiw[t] = - 1e-5 * g * (ρmeanw[1, 1, 1] - ρmeanw[1, 1, 69])

    ϕzc[t] = ρ₀ * sum(compute!(Field((wc - wbc / g) * vol))) 
    ϕzw[t] = ρ₀ * sum(compute!(Field((ww - wbw / g) * vol))) 
end

PEc = APEc + RPEc
PEw = APEw + RPEw

