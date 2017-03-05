/**
 * A worst-case O(N) time algorithm for the
 * total variation filter problem in one dimension.
 */

#include <vector>
#include <deque>

namespace tvf
{

/**
 * Given an input vector X of real numbers, the total variation
 * denoising of X is the vector Y minimising the function:
 *
 * \[ 1/2 * \| X - Y \|_2 + \lambda * \sum |Y_{i+1} - Y_{i}| \]
 *
 * This function calculates the total variation denoising of the
 * the input vector contained in the range [@p first, @last).
 * The values of the solution are written to the output iterator
 * @p output as values of type double. The function returns
 * an output iterator pointing to the end of the output range.
 *
 * Note: if any of the values [@p first, @p last) is NaN the
 * behaviour of the function is undefined.
 */
template<typename FwdIt, typename OutIt>
OutIt total_variation_denoise (
  FwdIt first,   FwdIt last,
  double lambda, OutIt output);

/*************************************************
 *
 * Inlined implementation
 *
 *************************************************/  
namespace detail
{
  
struct Point
{
  double x;
  double y;

  void operator+= (const Point &rhs)
  {
    x += rhs.x;
    y += rhs.y;
  }

  void operator-= (const Point &rhs)
  {
    x -= rhs.x;
    y -= rhs.y;
  }
};

inline double inner_product (const Point &u, const Point &v)
{
  return u.x * v.y - u.y * v.x;
}

/*
 * Add point on the right of convex/concave chain
 * removing all points not in the convex hull of the
 * new set.
 */
template<typename Q, typename Dot>
void convex_chain_extend (Q &chain, Point p, const Dot &dot)
{
  while (!chain.empty () && dot (p, chain.back ()) <= 0.)
  {
    p += chain.back ();
    chain.pop_back ();
  }

  chain.push_back (p);
}

/*
 * Symmetric
 */
template<typename InputIt, typename Dot>
InputIt
convex_chain_reduce (Point &p,
		     InputIt chain_begin,
		     const InputIt &chain_end,
		     const Dot &dot)
{
  while (chain_begin != chain_end &&
	 dot (p, *chain_begin) < 0)
  {
    p -= *chain_begin;
    ++chain_begin;
  }

  return chain_begin;
}

} // namespace detail {

template<typename InputIt, typename OutIt>
OutIt total_variation_denoise (InputIt first,
			       InputIt last,
			       const double lambda,
			       OutIt output)
{
  using namespace detail;

  /*
   * Geodesics to upper/lower bounding chains 
   */
  std::deque<Point> upper;
  std::deque<Point> lower;
  
  /*
   * Handle input of length <= 1 separately
   */
  if (first == last)
  {
    return output;
  }
  if (first + 1 == last)
  {
    *output++ = *first;
    return output;
  }
  else
  {
    const auto v = *first++;
    upper.push_back ({v + lambda, 1});
    lower.push_back ({v - lambda, 1});
  }
  
  /*
   * Logic to update upper/lower geodesics
   */
  const auto update = [&] (
    auto &u, auto &l, const Point &p, auto &&dot)
  {
    convex_chain_extend (u, p, dot);

    if (u.size () == 1 && !l.empty ())
    {
      const auto nf = convex_chain_reduce (
        u.front (), l.begin (), l.end (), dot);

      if (nf != l.begin ())
      {
	for (auto it = l.begin (); it != nf; ++it)
	{
	  output = std::fill_n (output, it->y, it->x / it->y);
	}

	l.erase (l.begin (), nf);
      }
    }
  };

  for (;;)
  {
    const auto value = *first++;
    const bool is_last = (first == last);

    const Point p {is_last ? (value - lambda) : value, 1};

    update (upper, lower, p, inner_product);
    
    if (is_last)
    {
      for (const auto &p : upper)
      {
	output = std::fill_n (output, p.y, p.x / p.y);
      }

      break;
    }
    else
    {
      update (lower, upper, p, [] (auto &&x, auto &&y)
      {
	return inner_product (y, x);
      });
    }
  }

  return output;
}

} // namespace tvf {
