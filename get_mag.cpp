#include "ipg_mag_vis.h"

double get_mag( std::vector<FFTResult> results ) {

    if (results.empty()) {
        std::cerr << "Error: No data found in 'fft_output'" << std::endl;
    }

    // Find the maximum magnitude and corresponding index
    FFTResult maxResult = results[0];
    for (const auto& result : results) {
		//~ std::cout << result.magnitude << std::endl;
        if (result.magnitude > maxResult.magnitude) {
            maxResult = result;
        }
    }

    return maxResult.magnitude;
}
