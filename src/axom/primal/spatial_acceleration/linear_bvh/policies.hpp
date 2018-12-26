#ifndef AXOM_PRIMAL_BVH_POLICIES_H_
#define AXOM_PRIMAL_BVH_POLICIES_H_

#include <RAJA/RAJA.hpp>

namespace axom
{
namespace primal
{
namespace bvh
{

#define BVH_CUDA_BLOCK_SIZE 128

#if defined(AXOM_USE_CUDA)
using raja_for_policy = RAJA::cuda_exec<BVH_CUDA_BLOCK_SIZE>;
using raja_reduce_policy = RAJA::cuda_reduce<BVH_CUDA_BLOCK_SIZE>;
using raja_atomic_policy = RAJA::atomic::cuda_atomic;
#elif defined(AXOM_USE_OPENMP)
using raja_for_policy = RAJA::omp_parallel_for_exec;
using raja_reduce_policy = RAJA::omp_reduce;
using raja_atomic_policy = RAJA::atomic::omp_atomic;
#else
using raja_for_policy = RAJA::seq_exec;
using raja_reduce_policy = RAJA::seq_reduce;
using raja_atomic_policy = RAJA::atomic::seq_atomic;
#endif
} /* namespace axom */
} /* namespace primal */
} /* namespace bvh */
#endif
