#include <random>                       /* for random number generator */

namespace axom
{
namespace mint
{

/*!
 * \brief Return a random double between min and max.
 * \param [in] min the lower bound.
 * \param [in] max the upper bound.
 */
inline double random_double( double min, double max )
{
  static std::random_device rd;
  static std::mt19937_64 mt( rd() );
  static std::uniform_real_distribution< double > dist(0.0, 1.0);
  
  double temp = dist(mt);
  return temp * ( max - min ) + min;
}

}	/* namespace mint */
}	/* namespace axom */