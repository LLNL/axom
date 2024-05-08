#ifndef SINA_ADIAK_HPP
#define SINA_ADIAK_HPP

/// @file

#include "axom/sina/config.hpp"
#ifdef AXOM_SINA_USE_ADIAK

#include <string>
#include <type_traits>

#include "axom/sina/include/ConduitUtil.hpp"
#include "axom/sina/include/Record.hpp"
#include "axom/sina/include/Run.hpp"

extern "C" {
#include "adiak_tool.h"
}

namespace axom
{
namespace sina
{

/**
 * The callback function to pass to Adiak in order to write collected data to a Sina Record.
 *
 * To register this with Adiak, you should call adiak_register_cb, passing in
 * a pointer to the record as the last value. Your code should look something
 * like this:
 *
 * \code
 *     Record record{ID{"my_id", IDType::Local}, "my_record_type"};
 *     adiak_register_cb(1, adiak_category_all, adiakSinaCallback, 0, &record);
 * \endcode
 *
 * Note that not everything that Sina can capture an be captured through the
 * Adiak API. For example, there is currently no support in Adiak to capture
 * anything like a CurveSet. As a result, to do that, you must hold on to
 * the Record object passed here as the opaque value and manipulate it directly.
 **/
void adiakSinaCallback(const char *name, adiak_category_t category, const char *subcategory, adiak_value_t *value, adiak_datatype_t *t, void *opaque_value);

}  // end sina namespace
}  // end axom namespace

#endif // AXOM_SINA_USE_ADIAK

#endif // SINA_ADIAK_HPP

