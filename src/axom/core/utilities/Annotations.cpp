#include "Annotations.hpp"

#include "axom/core/utilities/About.hpp"

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

#ifdef AXOM_USE_CALIPER
  #include "caliper/Caliper.h"
  #include "caliper/common/Attribute.h"
  #include "caliper/cali-manager.h"
  #ifdef AXOM_USE_MPI
    #include "caliper/cali-mpi.h"
  #endif
#endif

#ifdef AXOM_USE_CALIPER
namespace cali
{
extern Attribute region_attr;
}
#endif

namespace axom
{
namespace utilities
{
namespace annotations
{
static bool adiak_initialized = false;

#ifdef AXOM_USE_CALIPER
static cali::ConfigManager* cali_mgr {nullptr};
#endif

namespace detail
{
void initialize_adiak()
{
#ifdef AXOM_USE_ADIAK
  if(adiak_initialized)
  {
    return;
  }

  #ifdef AXOM_USE_MPI
  MPI_Comm communicator = MPI_COMM_WORLD;
  adiak::init((void*)&communicator);
  #else
  adiak::init(nullptr);
  #endif

  adiak::user();
  adiak::launchdate();
  adiak::executable();
  adiak::cmdline();
  adiak::clustername();
  adiak::jobsize();
  adiak::workdir();

  adiak::walltime();
  adiak::cputime();
  adiak::systime();

  adiak_initialized = true;
#endif
}

void initialize_caliper(const std::string& mode, int num_ranks)
{
#ifdef AXOM_USE_CALIPER
  bool multiprocessing = (num_ranks > 1);
  std::string configuration_service_list;
  cali::ConfigManager::argmap_t app_args;

  #ifdef AXOM_USE_MPI
  cali_mpi_init();
  #endif

  cali_mgr = new cali::ConfigManager();
  cali_mgr->add(mode.c_str(), app_args);

  for(const auto& kv : app_args)
  {
    const std::string& cali_mode = kv.first;
    const std::string& value = kv.second;

    if(cali_mode == "none")
    {
      // intentionally empty
    }
    else if(cali_mode == "report")
    {
      //report is an alias for the runtime-report Caliper configuration
      cali_mgr->add("runtime-report(output=stdout,calc.inclusive=true)");
    }
    else if(cali_mode == "counts")
    {
      configuration_service_list = "event:aggregate:";
      configuration_service_list += multiprocessing ? "mpireport" : "report";

      cali_config_preset("CALI_REPORT_CONFIG",
                         "SELECT count() "
                         "GROUP BY prop:nested "
                         "WHERE cali.event.end "
                         "FORMAT tree");

      cali_config_preset("CALI_MPIREPORT_CONFIG",
                         "SELECT   min(count) as \"Min count\", "
                         "         max(count) as \"Max count\", "
                         "         avg(count) as \"Avg count\", "
                         "         sum(count) as \"Total count\" "
                         "GROUP BY prop:nested "
                         "WHERE    cali.event.end "
                         "FORMAT   tree");
    }
    else if(cali_mode == "file")
    {
      configuration_service_list = "event:aggregate:timestamp:recorder";
    }
    else if(cali_mode == "trace")
    {
      configuration_service_list = "event:trace:timestamp:recorder";
    }
    else if(cali_mode == "gputx")
    {
  #if defined(AXOM_USE_CUDA)
      configuration_service_list = "nvtx";
  #elif defined(AXOM_USE_HIP)
      configuration_service_list = "roctx";
  #endif
    }
    else if(cali_mode == "nvtx" || cali_mode == "nvprof")
    {
      configuration_service_list = "nvtx";
    }
    else if(cali_mode == "roctx")
    {
      configuration_service_list = "roctx";
    }
    else if(!value.empty())
    {
      declare_metadata(cali_mode, value);
    }
  }

  if(!configuration_service_list.empty())
  {
    if(multiprocessing)
    {
      //This ensures mpi appears before any related service
      configuration_service_list = "mpi:" + configuration_service_list;
    }
    cali_config_preset("CALI_TIMER_SNAPSHOT_DURATION", "true");
    cali_config_preset("CALI_TIMER_INCLUSIVE_DURATION", "false");
    cali_config_preset("CALI_SERVICES_ENABLE",
                       configuration_service_list.c_str());
  }

  for(auto& channel : cali_mgr->get_all_channels())
  {
    channel->start();
  }

#endif  // AXOM_USE_CALIPER
}

}  // namespace detail

void initialize(const std::string& mode, int num_ranks)
{
  detail::initialize_adiak();
  detail::initialize_caliper(mode, num_ranks);

  declare_metadata("axom_version", axom::getVersion());
}

void finalize()
{
#ifdef AXOM_USE_ADIAK
  if(adiak_initialized)
  {
    adiak::fini();
  }
  adiak_initialized = false;
#endif
#ifdef AXOM_USE_CALIPER
  if(cali_mgr)
  {
    for(auto& channel : cali_mgr->get_all_channels())
    {
      channel->flush();
    }

    delete cali_mgr;
    cali_mgr = nullptr;
  }
#endif
}

void begin(const std::string& name)
{
#ifdef AXOM_USE_CALIPER
  cali::Caliper().begin(
    cali::region_attr,
    cali::Variant(CALI_TYPE_STRING, name.c_str(), name.length()));
#endif
}

void end(const std::string& /*name*/)
{
#ifdef AXOM_USE_CALIPER
  cali::Caliper().end(cali::region_attr);
#endif
}

}  // namespace annotations
}  // namespace utilities
}  // namespace axom