#include "dsptools.hpp"
#include <iostream>

// Brute force convolution
std::vector<double> dsp::conv(const std::vector<double>& x, const std::vector<double>& h, const char* mode)
{
	std::vector<double> y;
	
	// If the length of h is even, make it odd by adding an extra zero
	std::vector<double> h1;
	h1 = h;
	if (h1.size() % 2 == 0)
		h1.push_back(0);

	if(strcmp(mode, "same") == 0)
	{
		y.resize(x.size());

		for (int n = 0; n < x.size(); n++)
		{
			y[n] = 0;
			for (int m = -1*floor(h1.size()/2); m <= floor(h1.size()/2); m++)
			{
				double x_nm = 0;
				if(n-m >= 0 && n-m < x.size())
					x_nm = x[n-m];
				y[n] = y[n] + x_nm * h1[m+floor(h1.size()/2)];
			}
		}
	}
	else if(strcmp(mode, "full") == 0)
	{
		y.resize(x.size() + h1.size() - 1);
		
		for (int n = -1*floor(h1.size()/2); n < x.size()+floor(h1.size()/2); n++)
		{
			y[n+floor(h1.size()/2)] = 0;
			for (int m = -1*floor(h1.size()/2); m <= floor(h1.size()/2); m++)
			{
				double x_nm = 0;
				if(n-m >= 0 && n-m < x.size())
					x_nm = x[n-m];
				y[n+floor(h1.size()/2)] = y[n+floor(h1.size()/2)] + x_nm * h1[m+floor(h1.size()/2)];
			}
		}
	}
	else if(strcmp(mode, "valid") == 0)
	{
		y.resize(x.size() - h1.size() + 1);
		
		for (int n = floor(h1.size()/2); n < x.size()-floor(h1.size()/2); n++)
		{
			y[n-floor(h1.size()/2)] = 0;
			for (int m = -1*floor(h1.size()/2); m <= floor(h1.size()/2); m++)
			{
				double x_nm = 0;
				if(n-m >= 0 && n-m < x.size())
					x_nm = x[n-m];
				y[n-floor(h1.size()/2)] = y[n-floor(h1.size()/2)] + x_nm * h1[m+floor(h1.size()/2)];
			}
		}
	}
	
	return y;
}

// Converts the sampling rate of an input signal
std::vector<double> dsp::resample(const std::vector<double>& x, int Fs_old, int Fs_new)
{
	std::vector<double> y;

	if (Fs_old == Fs_new)
	{
		y = x;
		return x;
	}

	//----------- Step 1: Upsample the signal to the least common multiple frequency
	std::vector<double> x_up;
	double Fs_up = dsp::lcm(Fs_old, Fs_new);
	int up_factor = Fs_up / Fs_old;
	
	x_up.resize( (x.size()-1)*up_factor + 1 );
	for (int n = 0; n < x_up.size(); n++)
	{
		if (n % up_factor == 0)
			x_up[n] = x[n / up_factor];
		else
			x_up[n] = 0;
	}

	//----------- Step 2: Design a Low Pass Filter with cutoff fc, and apply it on x_up
	int down_factor = 1.0*Fs_up / Fs_new;
	double fc = 1.0 / (2.0*std::max(up_factor, down_factor)); // Normalized frequency (1 corresp. to Fs)
	int L = 20 * std::max(up_factor, down_factor);	// Filter order
	std::vector<double> lpf = dsp::design_fir_lowpass(L, fc);
	std::vector<double> x_lpf = dsp::conv(x_up, lpf, "same");
	for (int i = 0; i < x_lpf.size(); i++)
		x_lpf[i] = up_factor * x_lpf[i];

	//----------- Step 3: Downsample the signal x_lpf to get the output y
	y.resize( 1.0*(x.size()-1)*up_factor/down_factor + 1 );
	for (int n = 0; n < y.size(); n++)
	{
		y[n] = x_lpf[n * down_factor];
	}

	return y;
}

// Designs an FIR low pass filter using the window method.
// It uses a Hamming window.
std::vector<double> dsp::design_fir_lowpass(int L, double fc)
{
	// Make sure that L is even. If not, make it even.
	if (L % 2 != 0)
		L++;

	// Generate a truncated version of the ideal LPF impulse response
	std::vector<double> h;
	h.resize(L + 1);
	for (int i = 0; i < h.size(); i++)
	{
		if (i == L / 2)
			h[i] = 2 * fc;
		else
			h[i] = sin(2 * 3.1416 * fc * (i - L / 2)) / ((i - L / 2) * 3.1416);
	}

	// Generate a hamming window of length L+1
	std::vector<double> w;
	w.resize(L + 1);
	for (int i = 0; i < w.size(); i++)
		w[i] = 0.54 - 0.46*cos(2*3.1416*i/L);

	// Perform element-wise multiplication of h and w
	for (int i = 0; i < h.size(); i++)
		h[i] = h[i] * w[i];

	return h;
}

// Conventional filtering
std::vector<double> dsp::filter(std::vector<double> num, std::vector<double> den, const std::vector<double>& x)
{
	// Make both num and den coeffs same size
	if (num.size() > den.size())
		den.resize(num.size());
	else if (den.size() > num.size())
		num.resize(den.size());
	int N = num.size() - 1;	// Filter order

	// Make sure den[0] = 1
	if (den[0] != 1)
	{
		for (int i = N; i >= 0; i--)
		{
			num[i] = num[i] / den[0];
			den[i] = den[i] / den[0];
		}
	}

	// Initialize the state vector to zeros
	std::vector<double> s(N + 1, 0);

	// Filtering operation
	std::vector<double> y(x.size());
	for (int n = 0; n < x.size(); n++)
	{
		// Update s
		for (int i = N; i >= 1; i--)
			s[i] = s[i - 1];
		s[0] = x[n];
		for (int i = 1; i <= N; i++)
			s[0] = s[0] - den[i] * s[i];

		// Compute y[n]
		y[n] = 0;
		for (int i = 0; i <= N; i++)
			y[n] = y[n] + s[i] * num[i];
	}

	return y;
}

// FFT of a complex sequence using Cooley-Tukey algorithm
std::vector<std::complex<double> > dsp::fft(std::vector<std::complex<double> > x)
{
	int N = x.size();

	// Trivial base case
	if (x.size() == 1)
	{
		return(x);
	}

	// Make the size a power of 2 if it is already not
	int p = log2(N);
	if (pow(2, p) != N)
		x.resize(pow(2, p + 1), std::complex<double>(0, 0));
	N = x.size();

	// Reorder the elements
	std::vector<std::complex<double> > _x(N);
	for (int k = 0; k < N / 2; k++)
	{
		_x[k] = x[2 * k];
		_x[k + N / 2] = x[2 * k + 1];
	}
	x = _x;
	std::vector<std::complex<double> >().swap(_x);
	
	// Compute the FFTs of the two halves
	_fft_b(x.begin(), x.begin() + N/2);
	_fft_b(x.begin() + N/2, x.end());

	// Merge
	for (int k = 0; k <= N/2 - 1; k++)
	{
		std::complex<double> tmp = x[k];
		x[k] = tmp + (std::complex<double>(cos(2*3.1415*k/N), -1*sin(2*3.1415*k/N))) * x[k + N / 2];
		x[k + N / 2] = tmp - (std::complex<double>(cos(2 * 3.1415*k / N), -1 * sin(2 * 3.1415*k / N))) * x[k + N / 2];
	}

	return x;
}

// FFT of a real sequence (wrapper around complex FFT)
std::vector<std::complex<double> > dsp::fft(std::vector<double> x)
{
	std::vector<std::complex<double> > x_c(x.size());
	for (int i = 0; i < x.size(); i++)
		x_c[i] = std::complex<double>(x[i], 0);

	return dsp::fft(x_c);
}

// Internal function - in-place FFT
void dsp::_fft_b(std::vector<std::complex<double> >::iterator begin, std::vector<std::complex<double> >::iterator end)
{
	int N = std::distance(begin, end);

	// Trivial base case
	if (N == 1)
	{
		return;
	}

	// Reorder the elements
	std::vector<std::complex<double> > _x(N);
	for (int k = 0; k < N / 2; k++)
	{
		_x[k] = *(begin + 2 * k);
		_x[k + N / 2] = *(begin + 2 * k + 1);
	}
	for (int k = 0; k < N; k++)
		*(begin + k) = _x[k];
	std::vector<std::complex<double> >().swap(_x);

	// Compute the FFTs of the two halves
	_fft_b(begin, begin + N/2);
	_fft_b(begin + N/2, end);

	// Merge
	for (int k = 0; k <= N / 2 - 1; k++)
	{
		std::complex<double> tmp = *(begin + k);
		*(begin + k) = tmp + (std::complex<double>(cos(2 * 3.1416*k / N), -1 * sin(2 * 3.1416*k / N))) * (*(begin + k + N / 2));
		*(begin + k + N / 2) = tmp - (std::complex<double>(cos(2 * 3.1416*k / N), -1 * sin(2 * 3.1416*k / N))) * (*(begin + k + N / 2));
	}
}

// Finding peaks in signal
void dsp::find_peaks(const std::vector<double>& x, std::vector<double>& pks, std::vector<int>& lcs)
{
	int i = 1;
	while (i < x.size() - 1)
	{
		if (x[i] > x[i - 1] && x[i] > x[i + 1])		// Obviously a peak
		{
			pks.push_back(x[i]);
			lcs.push_back(i);
			i = i + 2;
		}
		else if (x[i] > x[i - 1] && x[i] == x[i + 1])	// Depends on slopes after i
		{
			// Keep going until a non-zero slope is reached
			int j;
			for (j = i + 1; j < x.size() - 1 && x[j] == x[j - 1]; j++);

			// If the slope is negative, it is a peak
			if (x[j] < x[j - 1])
			{
				pks.push_back(x[i]);
				lcs.push_back(i);
				i = j + 1;
			}
			else
				i = i + 1;
		}
		else
			i++;
	}
}

// General Math
int dsp::gcd(int a, int b)
{
	int c;
	while (a != 0) {
		c = a;
		a = b%a;
		b = c;
	}
	return b;
}

int dsp::lcm(int a, int b)
{
	return 1.0*(a*b)/dsp::gcd(a, b);
}
