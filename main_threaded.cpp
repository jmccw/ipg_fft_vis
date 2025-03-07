// comiling with 
//     g++ main_threaded.cpp fft_fftw.cpp get_mag.cpp -g -o main_threaded -lm -lfftw3


#include <iostream>
#include <cstdlib> // For std::atoi, rand()
#include <fstream>
#include <boost/array.hpp>
#include <vector>
#include <omp.h>

#include <boost/numeric/odeint.hpp>

#include "ipg_mag_vis.h"

using namespace std;
using namespace boost::numeric::odeint;

class System {
	
	public:
		// Define here the dimension of the system
		typedef boost::array< double , 3 > state_type;
		state_type x; 
		
		// Constants + Variables
		double alpha;
		double K;
		double P;
		double delta;
		double gamma_;
		
		std::vector<std::vector<double>> output;
		
		// Constructor | Reads in command line arguments, sets constants
		// This may need to be edited depending on the system.
		System(double alpha_in, double K_in, double delta_in, double gamma_in, double P_in, state_type x0)
			: alpha(alpha_in), K(K_in), delta(delta_in), gamma_(gamma_in), P(P_in), x(x0) {}
		
		
		
		void system_ode( const state_type &x , state_type &dxdt , double t )
		{
			dxdt[0] = x[2] * x[0] + K * std::cos(x[1]);
			dxdt[1] = -delta + alpha * x[2] - (this->K/x[0]) * std::sin(x[1]);
			dxdt[2] = gamma_ * (P - x[2] - (1 + 2*x[2]) * std::pow(x[0], 2));
		}

		void system_observer( const state_type &x , const double t ) 
		{
			this->output.push_back({t, x[0], x[1], x[2], x[0]*std::cos(x[1]), x[0]*std::sin(x[1])});
		}
		
		
		
		//~ void run() {  // Do not edit. Using lambdas to avoid scoping hassle.
			//~ integrate_const( boost::numeric::odeint::runge_kutta4<state_type>(),
				//~ [&](const state_type &x, state_type &dxdt, double t) {
					//~ system_ode(x, dxdt, t);
					//~ },
				//~ x, 0.0, 1200.0, 0.1,
				//~ [&](const state_type &x, double t) {
					//~ system_observer(x, t);
					//~ }
			//~ );
		//~ }
		
		void run() {  // Do not edit. Using lambdas to avoid scoping hassle.
			integrate( // dynamic stepsizing (really fast)
				[&](const state_type &x, state_type &dxdt, double t) {
					system_ode(x, dxdt, t);
					},
				x, 0.0, 1200.0, 0.1,
				[&](const state_type &x, double t) {
					system_observer(x, t);
					}
			);
		}

};

int main(int argc, char **argv)
{
	// Check if enough arguments are provided
    if (argc < 12) {
        std::cerr << "Usage: " << argv[0] << " <K_min> <K_max> <K_steps> <delta> <gamma_> <alpha> <R_0> <phi_0> <N_0> <P> <silent_mode>" << std::endl;
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
		
	// we need inputs for 
	// min&max of sweeping parameter
    
	double K_min = std::stod(argv[1]);
	double K_max = std::stod(argv[2]);
	int K_steps = std::stod(argv[3]);
	double delta_glob = std::stod(argv[4]);
	double gamma_glob = std::stod(argv[5]);
	double alpha_glob = std::stod(argv[6]);
	double R_0/*_upper*/ = std::stod(argv[7]);
	double phi_0/*_upper*/ = std::stod(argv[8]);
	double N_0/*_upper*/ = std::stod(argv[9]);
	double P_glob =  std::stod(argv[10]);
	int silent = (int)std::stod(argv[11]);
	
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
	//~ if(silent < 1){
		//~ std::cout << "K=" << K;
		//~ std::cout << " delta=" << delta;
		//~ std::cout << " gamma_=" << gamma_;
		//~ std::cout << " P=" << P;
		//~ std::cout << " alpha=" << alpha << std::endl;
	//~ }
	

	// data written to output vector    
    std::vector<std::vector<double>> max_magnitude;
    double delta_K = (K_max - K_min) / (double)(K_steps-1);
    
    #pragma omp parallel shared(max_magnitude) 
    {
		#pragma omp for 
		for(int j = 0; j < K_steps; j++){

			double K_thread = K_min + j*delta_K;
			System::state_type x0 = { 0.001 , 0.001, 0.001 };
			System simulation(alpha_glob, K_thread, delta_glob, gamma_glob, P_glob, x0);
			simulation.run();
				
			// Set fft data here - 10% of data is chopped to ignore settingly from initial conditions
			std::vector<double> fft_input;
			for(int i = (int)simulation.output.size()*0.75; i < (int)simulation.output.size(); i++) { 
				//for(int k = 0; k < (int)simulation.output[0].size(); k++) {
					fft_input.push_back(simulation.output[i][1]);
				//}
			}
			
			//std::vector<FFTResult> fft_result = fft_fftw( fft_input, silent ); // after integration do fft
			
			#pragma omp critical
			{
				// saving a alot of hassle here by just restricting to one fft at a time
				std::vector<FFTResult> fft_result = fft_fftw( fft_input, silent ); // after integration do fft
				max_magnitude.push_back( { K_thread, get_mag( fft_result ) } ); // returns largest mag value
				//~ std::cout<< magnitude <<std::endl;
			}
		}
	}
	
	//~ // write data out
	std::ofstream data_file("output.dat", std::ios::app); // Open output file
	for(int i = 0; i < max_magnitude.size(); i++){
		data_file << delta_glob << ' ';
		for(int j = 0; j < max_magnitude[0].size(); j++){
			data_file << max_magnitude[i][j] << ' ';
		}
		data_file << std::endl;
	}
	data_file.close(); 
	//~ // file will close automatically if it goes out of scope!
		
    if(silent < 1) std::cout << "success" << std::endl;
    return 0;
}
