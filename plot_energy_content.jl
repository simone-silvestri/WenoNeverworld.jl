using WenoNeverworld
using WenoNeverworld.Diagnostics
using GLMakie

prefix_simulation = "weno_fourth_ch"
dir = "/storage2/WenoNeverworldData/"
variables = ("b",)
stride = 5
fields = all_fieldtimeseries(prefix_simulation, dir; variables, checkpointer = true);

integrated_heat_content = Diagnostics.heat_content(fields[:b]; stride)

#kinetic_energy = Diagnostics.integral_kinetic_energy(fields[:u], fields[:v])
# transport = Diagnostics.ACC_transport(fields[:u])
#i = length(neverworld_fields[:u]) - 10
#KE = Diagnostics.KineticEnergy(neverworld_fields, i)
#vert_vorticity = Diagnostics.VerticalVorticity(neverworld_fields, i)
#rossby_radius also has one

APE = Diagnostics.integral_available_potential_energy(fields[:b]; stride)

#lines((APE))
GLMakie.lines(APE)
#GLMakie.lines(integrated_heat_content)
#GLMakie.lines(kinetic_energys)
