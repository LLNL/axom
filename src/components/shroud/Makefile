
tempdir := build/temp.linux-x86_64-2.7
testsdir := $(CURDIR)/tests

develop :
	python egg-setup.py develop

docs :
	python egg-setup.py build_sphinx

test :
	python egg-setup.py test
#	python -m unittest tests

compile : tutorial strings

tutorial strings :
	mkdir -p $(tempdir)/run-$@
	$(MAKE) \
	    -C $(tempdir)/run-$@ \
	    -f $(testsdir)/Makefile \
	     testsdir=$(testsdir) $@

test-compile : compile
	$(tempdir)/run-tutorial/tutorial
	$(tempdir)/run-strings/strings


.PHONY : develop docs test
.PHONY : compile test-compile tutorial strings
