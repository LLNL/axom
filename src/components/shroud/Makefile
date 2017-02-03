
top := $(CURDIR)
tempdir := build/temp.linux-x86_64-2.7
testsdir := $(top)/tests

include tests/defaults.mk

develop :
	python egg-setup.py develop

docs :
	python egg-setup.py build_sphinx

test :
	python egg-setup.py test
#	python -m unittest tests

# Pattern rule to make directories.
%/.. : ; $(at)test -d $(dir $@) || mkdir -p $(dir $@)

TESTDIRS = \
    $(tempdir)/run-tutorial/..\
    $(tempdir)/run-tutorial/python/.. \
    $(tempdir)/run-tutorial/lua/.. \
    $(tempdir)/run-strings/..

testdirs : $(TESTDIRS)

fortran : tutorial strings

tutorial strings : testdirs
	$(MAKE) \
	    -C $(tempdir)/run-$@ \
	    -f $(top)/tests/run-$@/Makefile \
	    top=$(top) $@

test-fortran : fortran
	$(tempdir)/run-tutorial/tutorial
	$(tempdir)/run-strings/strings

py-tutorial : testdirs
	$(MAKE) \
	    -C $(tempdir)/run-tutorial/python \
	    -f $(top)/tests/run-tutorial/python/Makefile \
	    top=$(top) all

test-python : py-tutorial
	export PYTHONPATH=$(top)/$(tempdir)/run-tutorial/python; \
	$(PYTHON_BIN) $(top)/tests/run-tutorial/python/test.py	


lua-tutorial : testdirs
	$(MAKE) \
	    -C $(tempdir)/run-tutorial/lua \
	    -f $(top)/tests/run-tutorial/lua/Makefile \
	    top=$(top) all

test-lua : lua-tutorial
#	export LUA_PATH=$(top)/$(tempdir)/run-tutorial/lua;
	cd $(top)/$(tempdir)/run-tutorial/lua; \
	$(LUA_BIN) $(top)/tests/run-tutorial/lua/test.lua

test-all : test-fortran test-python test-lua

.PHONY : develop docs test testdirs
.PHONY : fortran test-fortran tutorial strings
.PHONY : test-python py-tutorial
.PHONY : test-lua lua-tutorial
.PHONY : test-all
