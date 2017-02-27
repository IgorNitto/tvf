/**
 * A worst-case O(N) time algorithm for the
 * total variation filter problem in one dimension.
 */
#include <cmath>
#include <chrono>
#include <iostream>
#include <vector>
#include <deque>
#include <random>
#include <iomanip>

#include "tvf.hpp"

namespace
{

bool check (std::vector<double> in,
	    std::vector<double> filtered,
	    const double lambda)
{

  const auto float_equal = [] (double x, double y)
  {
    return (std::abs (x - y) < 1e-5) ||
           (std::abs (x - y) <
           1e-5 * ((std::abs (x) + std::abs (y)) / 2));
  };
  
  double s = 0;
  for (double &x : in)
  {
    x += s;
    s = x;
  }

  s = 0;
  for (double &x : filtered)
  {
    x += s;
    s = x;
  }

  std::cerr << std::setprecision(19);

  if (in.size () != filtered.size ())
  {
    std::cerr << "Inconsistent sizes" << std::endl;
    return false;
  }

  if (!float_equal (in.back (), filtered.back()))
  {
    std::cerr << "Mismatching overall sum " << in.back () << " " << filtered.back ();
    return false;
  }
  
  for (size_t j = 1; j < in.size () - 1; ++j)
  {
    const double m = (filtered[j - 1] + filtered[j + 1]) / 2;
    if (float_equal (filtered[j], m))
    {
      if (filtered[j] - in[j] > lambda * (1 + 1e-5) ||
	  filtered[j] - in[j] < -lambda * (1 + 1e-5))
      {
	std::cerr << "Vioated bounding tubes" <<
	          j << " " << filtered[j] << " " << in[j]
                  << " " << (in[j] - filtered[j]);

	return false;
      }
    }
    else if (filtered[j] < m)
    {
      if (!float_equal (filtered[j], in[j] + lambda))
      {
	std::cerr << j << " " << m << " " << filtered[j]
		  << " " << in[j] + lambda
		  << " " << in[j] - lambda << std::endl;
	
	return false;
      }
    }
    else if (!float_equal (filtered[j], in[j] - lambda))
    {
      std::cerr << j << " " << m << " " << filtered[j]
		<< " " << in[j] + lambda
		<< " " << in[j] - lambda << std::endl;
      
      return false;
    }
  }

  return true;
}

} // namespace {
  
int main (int, char**)
{
  std::default_random_engine source;
  std::uniform_real_distribution<double> gen (-0.1, 0.1);
  std::uniform_int_distribution<int> gen_int (-10, 10);

  std::vector<double> input;
  for (size_t j = 0; j < 100000; ++j)
  {
    const auto v = gen_int (source);
    for (size_t j = 0; j < 30; ++j)
    {
	input.push_back (v + gen (source));
    }
  }

  std::vector<double> output;
     
  tvf::total_variation_filter (
    input.begin (),
    input.end (),
    0.5,
    back_inserter (output));
   
  if (!check (input, output, 0.5))
  {
    std::cout << "Test failed\n";
    return EXIT_FAILURE;
  }
  else
  {
    std::cout << "Test passed\n";
  }
  return EXIT_SUCCESS;
}
