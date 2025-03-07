#ifndef IPG_MAG_VIS_H
#define IPG_MAG_VIS_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <complex>

struct FFTResult {
    size_t index;
    double real;
    double imag;
    double magnitude;
};

// Function declarations (not definitions)

std::vector<FFTResult> fft_fftw(std::vector<double> data, int silent);

double get_mag( std::vector<FFTResult> results );


#endif // IPG_MAG_VIS_H
