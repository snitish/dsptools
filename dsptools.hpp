#ifndef DSPTOOLS_HPP_INCLUDED
#define DSPTOOLS_HPP_INCLUDED

#include <vector>
#include <algorithm>
#include <complex>
#include <iterator>
#include <string.h>

//---------------------------------------------- Declaration of functions ----------------------------------------------
namespace dsp
{
	// Convolution
	std::vector<double> conv(const std::vector<double>& x, const std::vector<double>& h, const char* mode);
	
	// Resampling
	std::vector<double> resample(const std::vector<double>& x, int Fs_old, int Fs_new);

	// Filter design
	std::vector<double> design_fir_lowpass(int L, double fc);

	// Filtering
	std::vector<double> filter(std::vector<double> num, std::vector<double> den, const std::vector<double>& x);

	// FFT
	std::vector<std::complex<double> > fft(std::vector<double> x);
	std::vector<std::complex<double> > fft(std::vector<std::complex<double> > x);
	void _fft_b(std::vector<std::complex<double> >::iterator begin, std::vector<std::complex<double> >::iterator end);

	// Finding peaks in signal
	void find_peaks(const std::vector<double>& x, std::vector<double>& pks, std::vector<int>& lcs);

	// General math
	int gcd(int a, int b);
	int lcm(int a, int b);
}
//----------------------------------------------------------------------------------------------------------------------------


#endif