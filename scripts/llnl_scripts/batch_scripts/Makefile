# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
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
src-cz : src-cz-toss3 src-cz-blueos

src-cz-toss3 :
	ssh ruby      'cd $(here) && msub  msub_llnl_cz_toss3_src.sh'
	$(sleep)
src-cz-blueos :
	ssh lassen    'cd $(here) && bsub  < bsub_llnl_cz_blueos_src.sh'
	$(sleep)

##########
tpl-cz : sleep = sleep 2
tpl-cz : tpl-cz-toss3 tpl-cz-blueos
	$(info Results in /usr/WS1/axom/libs)

tpl-cz-toss3 :
	ssh ruby      'cd $(here) && msub  msub_llnl_cz_toss3_all_compilers.sh'
	$(sleep)
tpl-cz-blueos :
	ssh lassen    'cd $(here) && bsub  < bsub_llnl_cz_blueos_all_compilers.sh'
	$(sleep)

##########
src-rz : sleep = sleep 2
src-rz : src-rz-toss3 src-rz-blueos src-rz-blueosp9

src-rz-toss3 :
	ssh rzgenie   'cd $(here) && msub  msub_llnl_rz_toss3_src.sh'
	$(sleep)
src-rz-blueos :
	ssh rzansel   'cd $(here) && bsub  < bsub_llnl_rz_blueos_src.sh'
	$(sleep)
src-rz-blueosp9 :
	ssh rzmanta   'cd $(here) && bsub  < bsub_llnl_rz_blueos_src.sh'
	$(sleep)

##########
tpl-rz : sleep = sleep 2
tpl-rz : tpl-rz-toss3 tpl-rz-blueos tpl-rz-blueosp9
	$(info Results in /usr/WS1/axom/libs)

tpl-rz-toss3 :
	ssh rzgenie   'cd $(here) && msub  msub_llnl_rz_toss3_all_compilers.sh'
	$(sleep)
tpl-rz-blueos :
	ssh rzansel   'cd $(here) && bsub  < bsub_llnl_rz_blueos_all_compilers.sh'
	$(sleep)
tpl-rz-blueosp9 :
	ssh rzmanta   'cd $(here) && bsub  < bsub_llnl_rz_blueos_all_compilers.sh'
	$(sleep)

##########

# Remove batch logs
clean :
	rm -f *.out m.out.*


.PHONY : help
.PHONY : src-cz src-cz-toss3 src-cz-blueos
.PHONY : tpl-cz tpl-cz-toss3 tpl-cz-blueos
.PHONY : src-rz src-rz-toss3 src-rz-blueos src-rz-blueosp9
.PHONY : tpl-rz tpl-rz-toss3 tpl-rz-blueos tpl-rz-blueosp9
.PHONY : clean
