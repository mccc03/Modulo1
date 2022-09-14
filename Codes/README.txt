The file parameters.txt in the folder _data/input/ contains the input variables for the simulation, in the following order:

measures --> (int) number of measurements taken
resamplings --> (int) number of times the measurements are resampled
decorrel_len --> (int) number of times the matrix is updated before taking a measurement
Nlatt --> (int) length of the side of the matrix
init_flag --> (int) that initializes the matrix: if its value is (0) then all the spins are aligned (cold state), otherwise all the spins are random (hot state)
ext_field --> (double) value of the external magnetic field


The file ising_data_analysis.py stores the value of bta and Nlatt onto its output file, so keep this in mind when changing it during the simulations
