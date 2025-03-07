#include "ipg_mag_vis.h"

// Helper function to find the next power of 2
size_t nextPowerOf2(size_t n) {
    return std::pow(2, std::ceil(std::log2(n)));
}

//~ #define N 1024             // Number of samples
#define SAMPLE_RATE 1024.0 // Samples per second // this is here for completeness and to spoof algorithm, serves no purpose since we only care about magnitudes

void limit_precision(std::vector<double>& vec, int decimals) {
    double factor = std::pow(10.0, decimals);
    for (double& num : vec) {
        num = std::round(num * factor) / factor;
    }
}

double find_mean( std::vector<double> data ){
	double avg = 0.0;
	for( int i = 0; i < (int)data.size(); i++) {
		avg += data[i];
	}
	return avg / (double)data.size();
}

std::vector<FFTResult> fft_fftw(std::vector<double> data, int silent) {
	
	double avg = find_mean( data );
	for( int i = 0; i < (int)data.size(); i++) {
		data[i] -= avg;
	}
	limit_precision(data, 5); // limit precision to 5 decimal places so we're not capturing silly amounts of detail
    
    
    // Determine padded size
    int N = (int)data.size();
    int paddedSize = (int)nextPowerOf2(N);
    if(silent < 1){
		std::cout << "Original size: " << N << ", Padded size: " << paddedSize << std::endl;
	}

    // Pad data with zeros
    size_t begin_count = (paddedSize - N)/2;
    std::vector<double> paddedData(paddedSize, 0.0);
    for (size_t i = begin_count; i < N; ++i) {
        paddedData[i] = data[i];
    }
    int originalSize = N;
    
    // If I just don't do padding I dont get artefacts.
    // N = paddedSize; 
    // Keeping N at original size.

    if (data.empty()) {
        std::cerr << "FFT Error: No data found in input vector." << std::endl;
    }
    
    double *in;
    fftw_complex *out;
    fftw_plan plan;

    // Allocate memory
    in = (double*) fftw_malloc(sizeof(double) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2 + 1));
    plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);

    //~ double frequency = 50.0; // Hz
    //~ for(int i = 0; i < N; i++) {
        //~ in[i] = paddedData[i];
    //~ }
    
    for(int i = 0; i < originalSize; i++) {
        in[i] = data[i];
    }

    fftw_execute(plan);

	// return
	std::vector<FFTResult> fft_out;
    for(int i = 0; i < N/2 + 1; i++) {
		FFTResult result;
        result.real = out[i][0];
        result.imag = out[i][1];
        result.magnitude = sqrt(result.real*result.real + result.imag*result.imag)/(double)N;
        result.index = i * SAMPLE_RATE / N;
        fft_out.push_back(result);
        //~ printf("%f Hz: %f\n", current_freq, magnitude);
    }

    // Cleanup
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return fft_out;
}
