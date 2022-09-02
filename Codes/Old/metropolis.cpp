#include <iostream>
#include <fstream>
#include <cmath>

// Selecting computational cost
#define N 1024*1024

// Defining random number generator parameters
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

// Declaring simulation parameters
double x = 2.0;
double x0 = 5.0;
double sigma = 1.0;
double D = 0.1;
long int seed_start = 2;
long int * seed;

// Number of resamplings
unsigned int k = 256;

// Declaring miscellaneous variables
int i;
double mean, variance;
double position_array[N] = {0.0};
double * position;

// Functions

float ran2(long *idum);

double Gaussian(double x1);

double Mean(double * array);

void Metropolis();

void Bootstrap();

// Output files
std::ofstream out_file;
std::ofstream out_resample;

int main() {
  seed = &seed_start; // Assigning the seed
  position = &position_array[0];
  // Opening and overwriting the output file
  out_file.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/MetDistr.txt", std::ios::trunc);
  if (out_file.is_open()) {
    Metropolis();
    out_file.close();
  }

  out_resample.open("/home/exterior/Documents/Physics/MetodiNumerici/Modulo1/_data/mean_values.txt", std::ios::trunc);
  if(out_resample.is_open()){
    Bootstrap();
    out_resample.close();
  }

  return 0;
}

// Probability function
double Gaussian(double x1){
  return exp(-(x1-x0)*(x1-x0)/(2*sigma*sigma));
}

float ran2(long *idum){
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

// Metropolis algorithm
void Metropolis(){
  for(i=0;i<N;i++){
    double x_new = x+(2*ran2(seed)-1)*D;
    double u = ran2(seed);
    // Update check
    if(u<Gaussian(x_new)/Gaussian(x)){
      x = x_new;
    }
    *(position +i) = x;
    out_file << i << "\t" << x << "\n";
  }
  return;
}

//Bootstrap algorithm
void Bootstrap() {
  for(unsigned int p=0; p<N/k; p++){
    double sum = 0.0;
    int j = (int)(N*ran2(seed));
    for(unsigned int s=0; s<k; s++){
      if (j+k>N){
        sum = sum + *(position +(N-s-1));
        //*(position_resample +(p*l+s)) = *(position +(N-s-1));
      }
      else{
        sum = sum +*(position +(j+s));
        //*(position_resample +(p*l+s)) = *(position +(j+s));
      }
    }
    double mean_p = sum/(k);
    out_resample << p << "\t" << mean_p << "\n";
  }
  return;
}
