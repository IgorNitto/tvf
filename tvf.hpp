/**
 * A worst-case O(N) time algorithm for the
 * total variation filter problem in one dimension.
 */

#include <vector>
#include <deque>

namespace tvf
{

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

template<typename InputIt, typename Dot>
InputIt
convex_chain_reduce (Point &p,
		     InputIt chain_begin,
		     const InputIt &chain_end,
		     const Dot &dot)
{
  Point v {0., 0.};
  for (; chain_begin != chain_end; ++chain_begin)
  {
    auto next_v = v;
    next_v += *chain_begin;

    if (dot (p, next_v) > 0)
    {
      break;
    }

    v = next_v;
  }

  p -= v;
  return chain_begin;
}

} // namespace detail {

template<typename InputIt, typename OutIt>
OutIt total_variation_filter (InputIt first,
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

  if (std::distance (first, last) == 1)
  {
    *output++ = *first;
    return output;
  }
  else
  {
    const auto v = *first++;
    upper.push_back ({v + lambda/2, 1});
    lower.push_back ({v - lambda/2, 1});
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

    const Point p {is_last ? (value - lambda / 2) : value, 1};

    update (upper, lower, p, inner_product);

    if (is_last)
    {
      for (const auto &p : upper)
      {
	output = std::fill_n (output, p.y, p.x / p.y);
      }

      return output;
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
