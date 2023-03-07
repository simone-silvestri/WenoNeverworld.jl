# WenoNeverworld.jl

Modified neverworld2 setup from [(original neverworld2 setup)](https://egusphere.copernicus.org/preprints/2022/egusphere-2022-186/egusphere-2022-186.pdf) to 
1. study mesh convergence  
2. investigate numerics
3. calibrate and assess mesoscale parameterizations

# Domain and bathymetry

The starting point to simulate a _Neverworld_ setup is to construct a neverworld grid and add bathymetry. This is done with
```
grid = NeverworldGrid(arch, resolution, latitude = (-70, 70))
```

