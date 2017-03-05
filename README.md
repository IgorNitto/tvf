Total variation denoising
==========================================================

Solver for standard 1-D total variation denoising problem:

<img src="doc/TVD1.png">

The algorithm is exact up to floating point rounding error and run in
worst-case linear time in the input size. 

Reference: <https://en.wikipedia.org/wiki/Total_variation_denoising>

Usage
-----------------------------------

It is sufficient to include file `tvf.hpp`, no precompilation necessary.


```cpp

  #include "tvf.hpp"
  
  const std::vector<double> input {0, 10, 0};
  const double lambda = 0.5;
  
  std::vector<double> output (input.size ());

  tvf::total_variation_denoise (
    input.begin (), input.end (),
    lambda,
    output.begin ());  

```
