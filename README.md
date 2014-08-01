dsptools
========

DspTools is a collection of easy-to-use Digital Signal Processing algorithms implemented in C++. The aim is to make this library light, portable and easy to integrate with your code. All the functions in the library can work on std::vectors.

Usage
=====

DspTools is very simple to use. It's syntax is somewhat similar to that of MATLAB. Here is a sample program that demonstrates the basic usage:
```cpp
#include <iostream>
#include <cmath>
#include "dsptools.h"

int main(int argc, char** argv)
{
  // Generate a signal
  std::vector<double> x(100);
  for (int i=0; i<x.size(); i++)
    x[i] = 2*sin(2*3.1416*i/100) + cos(5*3.1415*i/100);
  
  // Finding the FFT
  std::vector<complex> x_fft = dsp::fft(x);
  
  // Filtering with an FIR filter
  std::vector<double> filter = {1, 1, 1, 1, 1};
  x = dsp::conv(x, filter, "same");
  
  return 0;
}

```
