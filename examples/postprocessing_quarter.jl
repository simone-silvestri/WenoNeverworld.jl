using WenoNeverworld
using WenoNeverworld.Diagnostics

all_fields = all_fieldtimeseries("files_four_centered/neverworld_quarter_centered_snapshots.jld2")
new_fields = limit_timeseries!(all_fields, all_fields[:b].times[end-20:end])

using WenoNeverworld.Diagnostics: average_spectra

spectraU = average_spectra(new_fields[:u], Colon(), 75:85; k = 65)