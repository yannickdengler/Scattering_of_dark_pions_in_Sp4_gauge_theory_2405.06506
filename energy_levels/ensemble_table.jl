using Pkg; Pkg.activate("./energy_levels/src_jl",io=devnull)
using I2julia

write_ensemble_list("output/logfiles.hdf5";outdir="output/tables",filename="pipi_fitintervals_default.csv")
