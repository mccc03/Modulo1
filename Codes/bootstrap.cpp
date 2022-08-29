#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <random>

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


///////////////
/* I/O files */
///////////////

ifstream input_Parameters;
ifstream energy_input;
ofstream energy_output;
ifstream magnetization_input;
ofstream magnetization_output;
ofstream susceptibility_output;
ofstream heatcapacity_output;
ifstream bta_input;


/////////////////////////////////////
/* Declaring resampling parameters */
/////////////////////////////////////

int measures;
int resamples;
int resampling_len;
double bta;
int Nlatt;


////////////////////////////////////
/* Initializing program variables */
////////////////////////////////////

double * energy_array;
double * magnetization_array;
double * energy_resampled_array;
double * magnetization_resampled_array;

/////////////////////////
/* Declaring functions */
/////////////////////////

float Ran2(long *idum); // Random number generator
double MeanValue(double * array, int len_array);
double Variance(double * array, int len_array);


//////////////////
/* Main program */
//////////////////

int main(){

    /* Initialize seed for RNG */
    seed = &seed_start;

    /* Read input parameters */
    input_Parameters.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/parameters.txt", ios::in);
    if(input_Parameters.is_open()){
        string line;
        string::size_type sz;
        getline(input_Parameters, line);
        measures = stoi(line,&sz);
        getline(input_Parameters, line);
        resampling_len = stoi(line,&sz);
        getline(input_Parameters, line);
        resamples = stoi(line,&sz);
        getline(input_Parameters, line);
        getline(input_Parameters, line);
        Nlatt = stoi(line,&sz);
        input_Parameters.close();
    }
    else{ // Error message
        cerr << "Could not read input parameters file";
    }


    bta_input.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/bta.txt", ios::in); // Opening beta.txt file, value stored in this file for easier access during recursion
    if(bta_input.is_open()){
        string line;
        string::size_type sz;
        getline(bta_input, line);
        bta = stod(line,&sz);
    }
    bta_input.close();

    /* Create dinamically allocated arrays - whatever that means - so that we don't run into segmentation faults for large values of measures variable */

    energy_array = new double[measures];
    energy_resampled_array = new double[measures];
    magnetization_array = new double[measures];
    magnetization_resampled_array = new double[measures];


    /* Open stored values for energy and magnetization */
    energy_input.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/energy.txt", ios::in);
    magnetization_input.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/magnetization.txt", ios::in);

    /* Read input files and store data onto arrays */
    if(energy_input.is_open() && magnetization_input.is_open()){
        for(int i=0;i<measures;i++){
            string line;
            string::size_type sz;
            getline(energy_input, line);
            *(energy_array+i) = stod(line,&sz);
            getline(magnetization_input, line);
            *(magnetization_array+i) = stod(line,&sz);
        }
    }
    else{ // Error message
        cerr << "Could not read input files";
    }

    /* Close input files */
    energy_input.close();
    magnetization_input.close();

    /* Open output files */
    energy_output.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/ene_res.txt", ios::trunc);
    magnetization_output.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/mag_res.txt", ios::trunc);
    susceptibility_output.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/sus_res.txt", ios::trunc);
    heatcapacity_output.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/heatc_res.txt", ios::trunc);

    /* Declare random int for resampling*/
    int j;

    /* Start bootstrap algorithm */
    if(energy_output.is_open() && magnetization_output.is_open() && susceptibility_output.is_open() && heatcapacity_output.is_open()){
        for(int k=0;k<resamples;k++){
            for(int p=0;p<measures/resampling_len;p++){
                for(int s=0;s<resampling_len;s++){
                    j = (int)(measures*Ran2(seed));
                    if(j+resampling_len>measures){
                        *(energy_resampled_array +p*resampling_len +s) = *(energy_array +(measures-s-1));
                        *(magnetization_resampled_array +p*resampling_len +s) = *(magnetization_array +(measures-s-1));
                    }
                    else{
                        *(energy_resampled_array +p*resampling_len +s) = *(energy_array +(j+s));
                        *(magnetization_resampled_array +p*resampling_len +s) = *(magnetization_array +(j+s));
                    }
                }
            }

            energy_output << k << "\t" << MeanValue(energy_resampled_array, measures) << "\n";
            magnetization_output << k << "\t" << MeanValue(magnetization_resampled_array, measures) << "\n";
            susceptibility_output << k << "\t" << bta*Nlatt*Nlatt*Variance(magnetization_resampled_array, measures) << "\n";
            heatcapacity_output << k << "\t" << bta*bta*Nlatt*Nlatt*Variance(energy_resampled_array, measures) << "\n";

        }
    }

    else{ // Error message
        cerr << "Could not open output files";
    }

    energy_output.close();
    magnetization_output.close();
    susceptibility_output.close();
    heatcapacity_output.close();

    /* Deallocating memory */
    delete[] energy_array;
    delete[] energy_resampled_array;
    delete[] magnetization_array;
    delete[] magnetization_resampled_array;


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


double MeanValue(double * array, int len_array){
  double sum = 0.0;
  for(int i=0;i<len_array;i++){
    sum = sum + array[i];
  }
  return ((double)sum)/len_array;
}


double Variance(double * array, int len_array){
    double sum_sq = 0.0;
    double mean = MeanValue(array,len_array);
    for(int i=0;i<len_array;i++){
        sum_sq = sum_sq + array[i]*array[i];
    }
    return (sum_sq/len_array) - mean*mean;
}
