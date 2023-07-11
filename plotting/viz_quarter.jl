using GLMakie, JDL2

# Load the data
@info "Loading data..."
path = pwd()
hfile = jldopen(path * "/weno_fourth.jld2", "r")
keys(hfile)
# pick out b field

η =  hfile["η"]["data"][:,:]

fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], xlabel = "x", ylabel = "y")
heatmap!(ax, η, colormap = :thermal)
display(fig)
save("weno_fourth_eta.png", fig)
