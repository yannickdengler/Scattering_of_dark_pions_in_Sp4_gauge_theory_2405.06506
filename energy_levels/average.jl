using Pkg; Pkg.activate("./energy_levels/src_jl")
using I2julia
using HDF5 

logfileshdf5 = "output/logfiles.hdf5"
write_averaged_hdf5_files(logfileshdf5)
