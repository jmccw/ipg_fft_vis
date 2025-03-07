// comiling with 
//     g++ main_threaded.cpp fft_fftw.cpp get_mag.cpp -g -o main_threaded -lm -lfftw3


#include <iostream>
#include <cstdlib> // For std::atoi, rand()
#include <fstream>
#include <boost/array.hpp>
#include <vector>

#include <boost/numeric/odeint.hpp>

#include "ipg_mag_vis.h"

using namespace std;
using namespace boost::numeric::odeint;

double alpha = 0.0;
double K = 0.195;
double P = 0.5;
double delta = 0.0;
double gamma_ = 0.05;

typedef boost::array< double , 3 > state_type;

std::vector<std::vector<double>> output;

void system_ode( const state_type &x , state_type &dxdt , double t )
{
    dxdt[0] = x[2] * x[0] + K * std::cos(x[1]);
    dxdt[1] = -delta + alpha * x[2] - (K/x[0]) * std::sin(x[1]);
    dxdt[2] = gamma_ * (P - x[2] - (1 + 2*x[2]) * std::pow(x[0], 2));
}

void system_observer( const state_type &x , const double t )
{
	output.push_back({t, x[0], x[1], x[2], x[0]*std::cos(x[1]), x[0]*std::sin(x[1])});
}

int main(int argc, char **argv)
{
	// Check if enough arguments are provided
    if (argc < 10) {
        std::cerr << "Usage: " << argv[0] << " <K> <delta> <gamma_> <alpha> <R_0> <phi_0> <N_0> <P> <silent_mode>" << std::endl;
        return 1;
    }
    
    // determining initial conditions (for general case bifs, random)
    //~ const long max_rand = 1000000L;
    //~ double lower_bound = 0;
    //~ double upper_bound = 10;
    //~ srandom(time(NULL));
 
	// Using random function from stack exchange:
    //~ double random_double = lower_bound
		//~ + (upper_bound - lower_bound)
		//~ * (random() % max_rand)
		//~ / max_rand;
    
	K = std::stod(argv[1]);
	delta = std::stod(argv[2]);
	gamma_ = std::stod(argv[3]);
	alpha = std::stod(argv[4]);
	double R_0/*_upper*/ = std::stod(argv[5]);
	double phi_0/*_upper*/ = std::stod(argv[6]);
	double N_0/*_upper*/ = std::stod(argv[7]);
	P =  std::stod(argv[8]);
	silent = (int)std::stod(argv[9]);
	
    //~ double R_0 = lower_bound + (R_0_upper - lower_bound)
                                 //~ * (random() % max_rand)
                                 //~ / max_rand;
    //~ double phi_0 = lower_bound + (phi_0_upper - lower_bound)
                                 //~ * (random() % max_rand)
                                 //~ / max_rand;
    //~ double N_0 = lower_bound + (N_0_upper - lower_bound)
                                 //~ * (random() % max_rand)
                                 //~ / max_rand;
                                 	
	//Parameter verification
	if(silent < 1){
		std::cout << "K=" << K;
		std::cout << " delta=" << delta;
		std::cout << " gamma_=" << gamma_;
		std::cout << " P=" << P;
		std::cout << " alpha=" << alpha << std::endl;
	}
	
	std::ofstream data_file_clear("output.dat", std::ios::trunc); // Clears the file
	data_file_clear.close(); // Close immediately after clearing
	// data written to output vector
    state_type x = { R_0 , phi_0, N_0 }; // initial conditions - here - need to add random seeding for initial conditions - more detail not getting picked up - read from cmd line input
    integrate( system_ode , x , 0.0 , 1200.0 , 0.01 , system_observer ); // boost function
        
    // Set fft data here - 10% of data is chopped to ignore settingly from initial conditions
    std::vector<double> fft_input;
    for(int i = (int)output.size()/10; i < (int)output.size(); i++) fft_input.push_back(output[i][1]);
    
    std::vector<FFTResult> fft_result = fft_fftw( fft_input ); // after integration do fft
    double max_magnitude = get_mag( fft_result ); // returns largest mag value
    //~ std::cout << max_magnitude << std::endl;
    // get magnitude of fft
    
    
    //~ std::ofstream data_file("output.dat", std::ios::app);
    //~ if (!data_file.is_open()) {
        //~ std::cerr << "Error: Unable to open file for writing." << std::endl;
    //~ }
    
    // write data out
    //~ for(int i = 0; i < output.size(); i++){
		//~ for(int j = 0; j < output[0].size(); j++){
			//~ data_file << output[i][j] << ' ';
		//~ }
		//~ data_file << std::endl;
	//~ }
	
    //~ data_file << t << ' ' << x[0] << ' ' << x[1] << ' ' << x[2] << std::endl;
    //~ data_file.close(); 
    // file will close automatically if it goes out of scope!
    
    if(silent < 1) std::cout << "success" << std::endl;
    return 0;
}
