loglistI2="./output/isospin_logfiles_I2_list"
loglist="./output/isospin_logfiles_list"

# activate correct julia environment and set up dependencies
julia -e 'using Pkg; Pkg.activate("./energy_levels/src_jl"); Pkg.instantiate()'
find ./input/ -name "out_scattering_I2" > $loglistI2
find ./input/ -name "out_spectrum" > $loglist
python3 energy_levels/HDF5.py $loglistI2

julia energy_levels/ensemble_table.jl
julia energy_levels/average.jl
julia energy_levels/write_fpi_correlator.jl
rm $loglistI2 $loglist

echo "Fit correlators: Obtain energy levels"
python3 energy_levels/fitting.py
julia energy_levels/energy_levels_table.jl

echo "Perform Lüscher analysis"
python3 scattering/scattering.py

echo Create Plots
python3 scattering/plotting.py
math -script scattering/sigma_v_I2.m
python3 scattering/plot_sigma.py
echo Done!

julia energy_levels/phaseshift_corrections.jl