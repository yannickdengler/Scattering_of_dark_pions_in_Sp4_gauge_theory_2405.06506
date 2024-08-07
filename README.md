# Scattering of dark pions in Sp(4) gauge theory
This repository contains the code used to prepare the plots and results included in [Scattering of dark pions in Sp(4) gauge theory [2405.06506]](https://arxiv.org/abs/2405.06506v1).

## Instructions: Running the analysis
- Install required dependencies (see below)
- Download the raw log files from the [Zenodo data release](10.5281/zenodo.12920978) (`isospin_logfiles`.zip), decompress the zip archive, and place it anywhere in the `input` directory
- Run the analysis using `bash main.sh` within the top-level directory
- The results used to prepare the manuscript [[2405.06506]](https://arxiv.org/abs/2405.06506v1) can then be found in `output/`

## Warning

The code in this repository has only been tested on the specific dataset provided here. It is not intended to be easily generalizable to arbitrary datasets. The analysis parameters are hard-coded in the directory `input/`.
The sections that are "commented out" should especially be treated with care.

## Requirements
- Mathematica 14
- Python 3.10 (see `requirements.txt` for the required packages)
- julia 1.10


### Installing python requirements using Conda

All python requirements may be installed from Conda by running

    conda env create -f environment.yml

On macOS on Apple silicon processors, it may be necessary to have Rosetta enabled, and to specify to use x86-64 packages, by running

    conda env create --platform osx-64 -f environment.yml

Before running the analysis,
the environment should be activated using

    conda activate Code_I2

## Logfile Parsing

The code "HDF5.py" in "Parsing" extracts the necessary information from the logfile obtained with the HiRep code and parses it to a HDF5 file. The code extracts the following information and saves it as:

- "logfile name": The name of the parsed logfile 
- "N_hits": The number of used sources
- "montecarlotimes": A list of the analysed trajectory numbers in Monte Carlo time.
- "filenames": A list of the file names of the analysed gauge configurations.
- "plaquette": The lattice averaged value of the plaquette for every gauge configuration
- "gauge_group": The gauge group
- "beta": The inverse bare coupling
- "m_1": The mass of the first fundamental fermion
- "m_2": The mass of the second fundamental fermion
- "N_L": The spatial extent of the lattice
- "N_T": The temporal extent of the lattice
- "operators": A list of the operators extracted from the logfile. Each operator has an entry for the real and the imaginary part.
- "correlators": A array that contains the correlation function of each measured operator. The array has four indices and has the size: ```num_Operators x num_soruces x num_Montecarlotimes x N_T```.

## Acknowledgments

YD and FZ have been supported the Austrian Science Fund research teams grant STRONG-DM (FG1). FZ has been supported by the STFC Grant No. ST/X000648/1. The computations have been performed on the Vienna Scientific Cluster (VSC4)

## External Data

The file "input/sigma_v_data_errors.csv" contains the data from FIG. 1 of [1]. It shows the velocity-weighted cross-section as a function of the mean velocity of dark matter in halos of dwarf-, low-surface brightness galaxies and galaxy clusters.

## References

[1] - Manoj Kaplinghat, Sean Tulin, and Hai-Bo Yu. Dark Matter Halos as Particle Colliders: Unified Solution to Small-Scale Structure Puzzles from Dwarfs to Clusters. [Phys. Rev. Lett., 116(4):041302, 2016. doi:10.1103/PhysRevLett.116.041302.](http://doi.org/10.1103/PhysRevLett.116.041302)
