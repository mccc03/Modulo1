#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;


/////////////////////////////
/* Defining RNG parameters */
/////////////////////////////

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

long int seed_start = 5;
long int * seed;


//////////////////////////
/* Defining spin matrix */
//////////////////////////

double ** spin_matrix;


/////////////////////////////////////
/* Declaring simulation parameters */
/////////////////////////////////////

int Nlatt; // Matrix dimension
/* Parameters to be read from parameters.txt */
int init_flag; // Starting state of the matrix, either hot (1), cold (0) or previous final state
int measures; // Number of measures
int decorrel_len; // Guess for decorrelation length
double ext_field; // External magnetic field's value
double bta; // Inverse of temperature

/* Pointers to the adjacent sites */
int * npp;
int * nmm;


///////////////
/* I/O files */
///////////////

ifstream input_Lattice;
ifstream input_Parameters;
ofstream output_Energy;
ofstream output_Magnetization;
ofstream output_Lattice;

ifstream bta_input;
ofstream bta_output;

ifstream old_Nlatt_input;
ofstream old_Nlatt_output;


/////////////////////////
/* Declaring functions */
/////////////////////////

float Ran2(long *idum); // Random number generator
void Geometry(int * movePlus, int * moveMinus); // Generates proximity arrays
void Lattice_init(double ** matrix, int flag, long int * seed); // Initializes matrix in a state defined by iflag
void Metropolis(double ** matrix, long int * seed); // Generates Markov chain
double Energy(double ** matrix); // Computes energy density
double Magnetization(double ** matrix); // Computes magnetization density
double Susceptibility(double ** matrix); // Computes susceptibility


//////////////////
/* Main program */
//////////////////

int main() {

    /* Reading the simulation parameters from parameters.txt */
    input_Parameters.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/input/parameters.txt", ios::in); // Opening file
    if (input_Parameters.is_open()) { // Check if file is open, then read

        string line;
        string::size_type sz;
        getline(input_Parameters, line);
        measures = stoi(line,&sz);
        getline(input_Parameters, line); // skip, number of resamplings
        getline(input_Parameters, line);
        decorrel_len = stoi(line,&sz);
        getline(input_Parameters, line);
        Nlatt = stoi(line,&sz);
        getline(input_Parameters, line);
        init_flag = stoi(line,&sz);
        getline(input_Parameters, line);
        ext_field = stod(line,&sz);

        input_Parameters.close();

    }

    else { // Error message
        cerr << "Unable to open parameters file.\n";
    }

    bta_input.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/input/bta.txt", ios::in); // Opening beta.txt file, value stored in this file for easier access during recursion
    if(bta_input.is_open()){
        string line;
        string::size_type sz;
        getline(bta_input, line);
        bta = stod(line,&sz);
    }
    bta_input.close();

    /* Initialize seed for RNG */
    seed = &seed_start;

    /* Create proximity arrays */
    int npp_array[Nlatt];
    int nmm_array[Nlatt];
    npp = npp_array;
    nmm = nmm_array;

    /* Compute proximity conditions */
    Geometry(npp,nmm);

    /* Allocating space for the spin matrix */
    spin_matrix = new double*[Nlatt];
    for(int count=0; count<Nlatt; count++){
        spin_matrix[count] = new double[Nlatt];
    }

    /* Initialize spin matrix */
    input_Lattice.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/lattice.txt", ios::in);
    old_Nlatt_input.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/old_Nlatt.txt", ios::in);
    if(input_Lattice.is_open() && old_Nlatt_input.is_open()){ // Check if file is open, then run
        Lattice_init(spin_matrix, init_flag, seed);
        input_Lattice.close();
    }
    else { // Error message
        cerr << "Unable to open Lattice file.\n";
    }

    /* Open output files */
    output_Energy.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/energy.txt", ios::trunc);
    output_Magnetization.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/magnetization.txt", ios::trunc);
    output_Lattice.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/lattice.txt", ios::trunc);

    /* Check if output files are open then run algorithm */
    if(output_Energy.is_open() && output_Magnetization.is_open()){

        /* Start Markov chain and take a measurement for each iteration */

        for(int l=0;l<measures;l++){

            /* Call Metropolis function decorrel_len times before taking a measurement*/
            for(int r=0;r<decorrel_len;r++){
                Metropolis(spin_matrix,seed);
            }

            /* Taking physical measurements */
            double ene = Energy(spin_matrix);
            double mag = Magnetization(spin_matrix);

            /* Writing measurements onto output files */
            output_Energy << ene << "\n";
            output_Magnetization << mag << "\n";

        }

    /* Closing output files */
    output_Energy.close();
    output_Magnetization.close();

    }

    else{// print error message
        cerr << "Unable to open output file(s).\n";
    }


    /* Prepare bta for next iteration */
    bta_output.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/input/bta.txt", ios::trunc);
    if(bta_output.is_open()){
        bta = bta+0.001;
        bta_output << bta << "\n";
    }
    bta_output.close();

    return 0;
}


///////////////
/* Functions */
///////////////

/* The Ran2 function is shamelessly taken from "Numerical recipes in C" */
float Ran2(long *idum){
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

/* This function creates two arrays that allow us to select the adjacent sites of a certain point, instead of computing them every time, which is more computationally demanding, while defining the "toroid" boundary conditions. */
void Geometry(int * movePlus, int * moveMinus){
    for(int count=0; count<Nlatt; count++){
        *(movePlus + count) = count+1;
        *(moveMinus + count) = count-1;
        *(movePlus + (Nlatt-1)) = 0; // Boundary
        *(moveMinus + 0) = Nlatt - 1; // Boundary
    }
  return;
}

/* The following function creates a matrix whose configuration depends on the value of init_flag */
void Lattice_init(double ** matrix, int flag, long int * seed){

    /* Initialize matrix whose elements are all 1.0 (cold) */
    if(flag==0){
        for(int row=0; row<Nlatt; row++){
            for(int column=0; column<Nlatt; column++){
                matrix[row][column] = 1.0;
            }
        }
    }

    /* Initialize matrix whose elements are random (hot) */
    else if(flag==1){
        for(int row=0; row<Nlatt; row++){
            for(int column=0; column<Nlatt; column++){
                float random_number = Ran2(seed);
                if(random_number<0.5){
                    matrix[row][column] = -1.0;
                }
                else{
                    matrix[row][column] = 1.0;
                }
            }
        }
    }

    /* Initialize matrix whose elements are read from the last iteration
    of the algorithm */
    else{
        string line;
        string::size_type sz;
        getline(old_Nlatt_input, line);
        int old_Nlatt = stod(line,&sz);
        if(old_Nlatt==Nlatt){
            for(int row=0; row<Nlatt; row++){
                for(int column=0; column<Nlatt; column++){
                    string line;
                    string::size_type sz;
                    getline(input_Lattice, line);
                    matrix[row][column] = stod(line,&sz);
                }
            }
        }
        else{
            cerr << "Can't load matrix of a different size.\n";
        }
    }

    return;
}

/* The Metropolis function defines the Markov chain. */
void Metropolis(double ** matrix, long int * seed){

    for(int ivol=0; ivol<(Nlatt*Nlatt); ivol++){

        // Selecting a random site of the matrix
        int i = (int)(round(Ran2(seed)*(Nlatt-1)));
        int j = (int)(round(Ran2(seed)*(Nlatt-1)));

        // Using the results of Geometry to select the adjacent sites
        int ip = *(npp + i);
        int im = *(nmm + i);
        int jp = *(npp +j);
        int jm = *(nmm +j);

        // This is where "physics" kicks in, we need to compute the "Hamiltonian"
        // of the system and update the matrix based on that value
        double force = matrix[i][jp]+matrix[i][jm]+matrix[ip][j]+matrix[im][j];
        force = bta*(force+ext_field);

        double phi = matrix[i][j]; // Current spin value

        // Probability ratio with inverted spin
        double prob_ratio = exp(-2.0*phi*force);

        // Last thing we need is the Metropolis test
        float u = Ran2(seed);
        if(u<prob_ratio){
            matrix[i][j] = -phi;
        }
    }

  return;
}

/* The Energy function explores the whole matrix and computes the energy density */
double Energy(double ** matrix){

    // This will be inside a cycle, so I have to reset energy_c
    double energy_c = 0.0;

    // Cycle over the entire matrix
    for(int i=0; i<Nlatt; i++){
        for(int j=0; j<Nlatt; j++){
            int ip = *(npp + i);
            int im = *(nmm + i);
            int jp = *(npp +j);
            int jm = *(nmm +j);
            double force = matrix[i][jp]+matrix[i][jm]+matrix[ip][j]+matrix[im][j];
            /* 0.5 factor is needed because I am counting
            each site twice during summation*/
            energy_c = energy_c - 0.5*force*matrix[i][j];
            energy_c = energy_c - ext_field*matrix[i][j];
        }
    }

    // What I need is actually the energy density
    return energy_c/(double)(Nlatt*Nlatt);
}

/* The Magnetization function explores the matrix and returns the abs value of its normalized element-wise sum */
double Magnetization(double ** matrix){

    // This will be inside a cycle, so I have to reset magnetization_c
    double magnetization_c = 0.0;

    // Cycle over the entire matrix
    for(int i=0; i<Nlatt; i++){
        for(int j=0; j<Nlatt; j++){
            magnetization_c = magnetization_c + matrix[i][j];
        }
    }
    // I actually need the average magnetization
    return abs(magnetization_c)/(double)(Nlatt*Nlatt);
}


