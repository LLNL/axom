// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/Annotations.hpp"

#include "axom/fmt.hpp"

#ifdef AXOM_USE_ADIAK
  #include <adiak_tool.h>
#endif

namespace
{
#ifdef AXOM_USE_ADIAK
static std::string adiak_value_as_string(adiak_value_t *val, adiak_datatype_t *t)
{
  // Implementation adapted from adiak user docs

  if(!t)
  {
    return "ERROR";
  }

  auto get_vals_array = [](adiak_datatype_t *t, adiak_value_t *val, int count) {
    std::vector<std::string> s;
    for(int i = 0; i < count; i++)
    {
      adiak_value_t subval;
      adiak_datatype_t *subtype;
      adiak_get_subval(t, val, i, &subtype, &subval);
      s.push_back(adiak_value_as_string(&subval, subtype));
    }
    return s;
  };

  switch(t->dtype)
  {
  case adiak_type_unset:
    return "UNSET";
  case adiak_long:
    return axom::fmt::format("{}", val->v_long);
  case adiak_ulong:
    return axom::fmt::format("{}", static_cast<unsigned long>(val->v_long));
  case adiak_longlong:
    return axom::fmt::format("{}", val->v_longlong);
  case adiak_ulonglong:
    return axom::fmt::format("{}",
                             static_cast<unsigned long long>(val->v_longlong));
  case adiak_int:
    return axom::fmt::format("{}", val->v_int);
  case adiak_uint:
    return axom::fmt::format("{}", static_cast<unsigned int>(val->v_int));
  case adiak_double:
    return axom::fmt::format("{}", val->v_double);
  case adiak_date:
    // holds time in seconds since epoch
    return axom::fmt::format(
      "{:%a %d %b %Y %T %z}",
      std::chrono::system_clock::time_point {std::chrono::seconds {val->v_long}});
  case adiak_timeval:
  {
    const auto *tv = static_cast<struct timeval *>(val->v_ptr);
    return axom::fmt::format("{:%S} seconds:timeval",
                             std::chrono::seconds {tv->tv_sec} +
                               std::chrono::microseconds {tv->tv_usec});
  }
  case adiak_version:
    return axom::fmt::format("{}:version", static_cast<char *>(val->v_ptr));
  case adiak_string:
    return axom::fmt::format("{}", static_cast<char *>(val->v_ptr));
  case adiak_catstring:
    return axom::fmt::format("{}:catstring", static_cast<char *>(val->v_ptr));
  case adiak_path:
    return axom::fmt::format("{}:path", static_cast<char *>(val->v_ptr));
  case adiak_range:
    return axom::fmt::format("{}",
                             axom::fmt::join(get_vals_array(t, val, 2), " - "));
  case adiak_set:
    return axom::fmt::format(
      "[{}]",
      axom::fmt::join(get_vals_array(t, val, adiak_num_subvals(t)), ", "));
  case adiak_list:
    return axom::fmt::format(
      "{{{}}}",
      axom::fmt::join(get_vals_array(t, val, adiak_num_subvals(t)), ", "));
  case adiak_tuple:
    return axom::fmt::format(
      "({})",
      axom::fmt::join(get_vals_array(t, val, adiak_num_subvals(t)), ", "));
  }
}

static void get_namevals_as_map(const char *name,
                                int AXOM_UNUSED_PARAM(category),
                                const char *AXOM_UNUSED_PARAM(subcategory),
                                adiak_value_t *value,
                                adiak_datatype_t *t,
                                void *opaque_value)
{
  // add name/value to adiak metadata map
  using kv_map = std::map<std::string, std::string>;
  auto &metadata = *static_cast<kv_map *>(opaque_value);
  metadata[name] = adiak_value_as_string(value, t);
}

#endif  // AXOM_USE_ADIAK
}  // namespace

TEST(utils_annotations, initialize_finalize)
{
  axom::utilities::annotations::initialize();

  axom::utilities::annotations::finalize();

  SUCCEED();
}

TEST(utils_annotations, print_adiak_metadata)
{
  axom::utilities::annotations::initialize();

#ifdef AXOM_USE_ADIAK
  std::map<std::string, std::string> metadata;
  adiak_list_namevals(1, adiak_category_all, get_namevals_as_map, &metadata);

  std::cout << "Adiak metadata: \n";
  for(const auto &kv : metadata)
  {
    std::cout << axom::fmt::format("- {}: {}\n", kv.first, kv.second);
  }
#endif

  axom::utilities::annotations::finalize();

  SUCCEED();
}
