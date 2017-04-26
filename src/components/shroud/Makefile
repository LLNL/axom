
top := $(CURDIR)

PYTHON := $(shell which python)
# 2.7
PYTHON_VER := $(shell $(PYTHON) -c "import sys;sys.stdout.write('{version[0]}.{version[1]}'.format(version=sys.version_info))")
# linux-x86_64
PLATFORM := $(shell $(PYTHON) -c "import sys, sysconfig;sys.stdout.write(sysconfig.get_platform())")

LUA = $(shell which lua)

# build/temp-linux-x86_64-3.6
tempdir := build/temp.$(PLATFORM)-$(PYTHON_VER)
testsdir := $(top)/tests

include $(top)/tests/defaults.mk

develop :
	python setup.py develop

docs :
	python setup.py build_sphinx

test :
	python setup.py test
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
	    PYTHON=$(PYTHON) top=$(top) all

test-python : py-tutorial
	export PYTHONPATH=$(top)/$(tempdir)/run-tutorial/python; \
	$(PYTHON_BIN) $(top)/tests/run-tutorial/python/test.py	

lua-tutorial : testdirs
	$(MAKE) \
	    -C $(tempdir)/run-tutorial/lua \
	    -f $(top)/tests/run-tutorial/lua/Makefile \
	    LUA=$(LUA) top=$(top) all

test-lua : lua-tutorial
#	export LUA_PATH=$(top)/$(tempdir)/run-tutorial/lua;
	cd $(top)/$(tempdir)/run-tutorial/lua; \
	$(LUA_BIN) $(top)/tests/run-tutorial/lua/test.lua

test-all : test-fortran test-python test-lua

test-clean :
	rm -rf $(tempdir)

print-debug:
	@echo LUA=$(LUA)
	@echo PYTHON=$(PYTHON)
	@echo PYTHON_PREFIX=$(PYTHON_PREFIX)
	@echo PYTHON_VER=$(PYTHON_VER)
	@echo PLATFORM=$(PLATFORM)
	@echo tempdir=$(tempdir)

.PHONY : develop docs test testdirs
.PHONY : fortran test-fortran tutorial strings
.PHONY : test-python py-tutorial
.PHONY : test-lua lua-tutorial
.PHONY : test-all test-clean
.PHONY : print-debug
