# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#
# Run batch jobs to build third-party-libraires or source for
# multiple configurations and hosts
#
# Since directories are named by time (HH_MM_SS), add a sleep between
# submitting jobs to minimize getting duplicate directory names.
# However, names may still be duplicated depending on when the batch
# system actually starts the job.

here := $(CURDIR)

help :
	$(info Start batch jobs on a machines of each SYS_TYPE from a single machine)
	$(info Must be run once from a RZ machine and once from a CZ machine)
	$(info since CURDIR will not be the same.)
	$(info )
	$(info Build src on CZ or RZ)
	$(info    make src-cz)
	$(info    make src-rz)
	$(info Build TPL and src on CZ or RZ)
	$(info    make tpl-cz)
	$(info    make tpl-rz)
	$(info )

##########
src-cz : sleep = sleep 2
src-cz : src-cz-toss3 src-cz-blueos src-cz-bgq

src-cz-toss3 :
	ssh quartz    'cd $(here) && msub -d $(here) msub_llnl_cz_toss3_src.sh'
	$(sleep)
src-cz-blueos :
	ssh ray       'cd $(here) && bsub  < bsub_llnl_cz_blueos_src.sh'
	$(sleep)
src-cz-bgq :
	ssh vulcan    'cd $(here) && msub -d $(here) msub_llnl_cz_bgq_src.sh'
	$(sleep)

##########
tpl-cz : sleep = sleep 2
tpl-cz : tpl-cz-toss3 tpl-cz-blueos tpl-cz-bgq
	$(info Results in /usr/WS1/axom/thirdparty_libs/builds)

tpl-cz-toss3 :
	ssh quartz    'cd $(here) && msub -d $(here) msub_llnl_cz_toss3_all_compilers.sh'
	$(sleep)
tpl-cz-blueos :
	ssh ray       'cd $(here) && bsub  < bsub_llnl_cz_blueos_all_compilers.sh'
	$(sleep)
tpl-cz-bgq :
	ssh vulcan    'cd $(here) && msub -d $(here) msub_llnl_cz_bgq_all_compilers.sh'
	$(sleep)

##########
src-rz : sleep = sleep 2
src-rz : src-rz-chaos5 src-rz-toss3 src-rz-blueos src-rz-bgq

src-rz-chaos5 :
	ssh rzalastor 'cd $(here) && sbatch sbatch_llnl_rz_chaos5_src.sh'
	$(sleep)
src-rz-toss3 :
	ssh rzgenie   'cd $(here) && msub -d $(here) msub_llnl_rz_toss3_src.sh'
	$(sleep)
src-rz-blueos :
	ssh rzmanta   'cd $(here) && bsub  < bsub_llnl_rz_blueos_src.sh'
	$(sleep)
src-rz-bgq :
	ssh rzuseq    'cd $(here) && msub -d $(here) msub_llnl_rz_bgq_src.sh'
	$(sleep)

##########
tpl-rz : sleep = sleep 2
tpl-rz : tpl-rz-chaos5 tpl-rz-toss3 tpl-rz-blueos tpl-rz-bgq
	$(info Results in /usr/WS1/axom/thirdparty_libs/builds)

tpl-rz-chaos5 :
	ssh rzalastor 'cd $(here) && sbatch sbatch_llnl_rz_chaos5_all_compilers.sh'
	$(sleep)
tpl-rz-toss3 :
	ssh rzgenie   'cd $(here) && msub -d $(here) msub_llnl_rz_toss3_all_compilers.sh'
	$(sleep)
tpl-rz-blueos :
	ssh rzmanta   'cd $(here) && bsub  < bsub_llnl_rz_blueos_all_compilers.sh'
	$(sleep)
tpl-rz-bgq :
	ssh rzuseq    'cd $(here) && msub -d $(here) msub_llnl_rz_bgq_all_compilers.sh'
	$(sleep)

##########

# Remove batch logs
clean :
	rm -f *.out m.out.*


.PHONY : help
.PHONY : src-cz src-cz-toss3 src-cz-blueos src-cz-bgq
.PHONY : tpl-cz tpl-cz-toss3 tpl-cz-blueos tpl-cz-bgq
.PHONY : src-rz src-rz-chaos5 src-rz-toss3 src-rz-blueos src-rz-bgq
.PHONY : tpl-rz tpl-rz-chaos5 tpl-rz-toss3 tpl-rz-blueos tpl-rz-bgq
.PHONY : clean
