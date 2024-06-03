# Scattering of dark pions in Sp(4) gauge theory
This repository contains the code used to prepare the plots and results included in [Scattering of dark pions in Sp(4) gauge theory [2405.06506]](https://arxiv.org/abs/2405.06506v1).

## Instructions: Running the analysis
- Install required dependencies (see below)
- Download the raw log files from the [Zenodo data release]() and place it in anywhere in the `input` directory
- Run the analysis using `bash main.sh` within the top-level directory
- The results used to prepare the manuscript [[2405.06506]](https://arxiv.org/abs/2405.06506v1) can then be found in `output/`

## Warning

The code in this repository has only been tested on the specific dataset provided here. It is not intended to be easily generalizable to arbitrary datasets. The analysis parameters are hard-coded in the directory `input/`.

## Requirements
- Python 3.8 (see `requirements.txt` for the required packages)
- julia 1.9

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
