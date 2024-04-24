#ifndef SINA_ADIAK_HPP
#define SINA_ADIAK_HPP

/// @file

#include "sina/config.hpp"
#if SINA_BUILD_ADIAK_BINDINGS

#include <string>
#include <type_traits>

#include "sina/ConduitUtil.hpp"
#include "sina/Record.hpp"
#include "sina/Run.hpp"

extern "C" {
#include "adiak_tool.h"
}

namespace sina {

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

}

#endif // SINA_BUILD_ADIAK_BINDINGS

#endif // SINA_ADIAK_HPP

