The file parameters.txt in the folder _data/input/ contains the input variables for the simulation, in the following order:

measures --> (int) number of measurements taken
resamplings --> (int) number of times the measurements are resampled
decorrel_len --> (int) number of times the matrix is updated before taking a measurement
Nlatt --> (int) length of the side of the matrix
init_flag --> (int) that initializes the matrix: if its value is (0) then all the spins are aligned (cold state), if its value is (1) all the spins are random (hot state), otherwise the code loads the last lattice configuration of the preious simulation
ext_field --> (double) value of the external magnetic field
thermal_len --> (int) measures to be discarded

At the end of ising.cpp, bta is increased by 0.005

The Shell file start_simulation.sh is used to iterate the simulation for different bta's
