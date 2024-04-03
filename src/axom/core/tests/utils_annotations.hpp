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
  #include "adiak_tool.h"
#endif

#ifdef AXOM_USE_CALIPER
  #include "caliper/cali.h"
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
  axom::utilities::annotations::initialize("none", 1);

  axom::utilities::annotations::finalize();

  SUCCEED();
}

TEST(utils_annotations, print_adiak_metadata)
{
  using MetadataMap = std::map<std::string, std::string>;

  axom::utilities::annotations::initialize("none", 1);

#ifdef AXOM_USE_ADIAK
  {
    MetadataMap default_metadata;
    adiak_list_namevals(1,
                        adiak_category_all,
                        get_namevals_as_map,
                        &default_metadata);

    // some default metadata is present, e.g. jobsize
    EXPECT_TRUE(default_metadata.find("user") != default_metadata.end());

    // but other metadata is not present before we register it
    EXPECT_TRUE(default_metadata.find("custom_int_metadata") ==
                default_metadata.end());
    EXPECT_TRUE(default_metadata.find("custom_str_metadata") ==
                default_metadata.end());
  }

  // register some metadata
  AXOM_ANNOTATE_METADATA("custom_str_metadata", "arbitrary string", "");
  AXOM_ANNOTATE_METADATA("custom_int_metadata", 42, "");

  // check that registered metadata is now present
  MetadataMap updated_metadata;
  adiak_list_namevals(1, adiak_category_all, get_namevals_as_map, &updated_metadata);

  // some default metadata is present
  EXPECT_TRUE(updated_metadata.find("user") != updated_metadata.end());

  // but other metadata is not present before we register it
  EXPECT_TRUE(updated_metadata.find("custom_int_metadata") !=
              updated_metadata.end());
  EXPECT_TRUE(updated_metadata.find("custom_str_metadata") !=
              updated_metadata.end());

  EXPECT_EQ(std::to_string(42), updated_metadata["custom_int_metadata"]);
  EXPECT_EQ("arbitrary string", updated_metadata["custom_str_metadata"]);

  // print the key-value pairs
  std::cout << "Adiak metadata: \n";
  for(const auto &kv : updated_metadata)
  {
    std::cout << axom::fmt::format("- {}: {}\n", kv.first, kv.second);
  }
#endif

  axom::utilities::annotations::finalize();
}

TEST(utils_annotations, check_modes)
{
  EXPECT_TRUE(axom::utilities::annotations::detail::check_mode("none"));

  EXPECT_TRUE(axom::utilities::annotations::detail::check_mode("counts"));
  EXPECT_TRUE(axom::utilities::annotations::detail::check_mode("file"));
  EXPECT_TRUE(axom::utilities::annotations::detail::check_mode("trace"));
  EXPECT_TRUE(axom::utilities::annotations::detail::check_mode("report"));

  EXPECT_TRUE(axom::utilities::annotations::detail::check_mode("gputx"));
  EXPECT_TRUE(axom::utilities::annotations::detail::check_mode("nvprof"));
  EXPECT_TRUE(axom::utilities::annotations::detail::check_mode("nvtx"));
  EXPECT_TRUE(axom::utilities::annotations::detail::check_mode("roctx"));

  EXPECT_FALSE(
    axom::utilities::annotations::detail::check_mode("_unregistered_"));
}

TEST(utils_annotations, modes)
{
  // TODO: Add tests for the following modes: "none", "report", "counts", "file", "trace"
  //       Must be done is separate runs since caliper does not cleanly tear down

  const std::string mode = "counts";
  std::cout << "Testing caliper service '" << mode << "'\n";

  axom::utilities::annotations::initialize(mode, 1);

#ifdef AXOM_USE_CALIPER
  {
    AXOM_ANNOTATE_BEGIN("my region");

    for(int i = 0; i < 10; ++i)
    {
      AXOM_ANNOTATE_SCOPE("inner1");
    }

    for(int i = 0; i < 100; ++i)
    {
      AXOM_ANNOTATE_SCOPE("inner2");
    }
    AXOM_ANNOTATE_END("my region");
  }
#endif

  axom::utilities::annotations::finalize();

  SUCCEED();
}

TEST(utils_annotations, print_help)
{
  std::cout << "Caliper help string: \n"
            << axom::utilities::annotations::detail::help_string() << std::endl;

  SUCCEED();
}
