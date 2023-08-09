using WenoNeverworld
using WenoNeverworld.Diagnostics
using GLMakie

prefix_simulation = "weno_one_ch"
dir = "/storage2/WenoNeverworldData/"
variables = ("b",)
stride = 10
fields = all_fieldtimeseries(prefix_simulation, dir; variables, checkpointer = true);
heat_content = Diagnostics.heat_content(fields[:b]; stride)
# kinetic_energy = Diagnostics.integral_kinetic_energy(fields[:u], fields[:v])
# transport = Diagnostics.ACC_transport(fields[:u])

lines(sort(heat_content))