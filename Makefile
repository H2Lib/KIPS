
# ------------------------------------------------------------
# Components of the main library
# ------------------------------------------------------------

KIPS_CORE0 = \
	Library/basic.c \
	Library/settings.c \
	Library/parameters.c \
	Library/interpolation.c \
  	Library/blas.c 

KIPS_CORE1 = \
	Library/avector.c \
	Library/rvector.c \
	Library/amatrix.c 

KIPS_CORE2 = \
	Library/quaternion.c \
	Library/spatialcluster.c \
	Library/spatialgeometry.c \
	Library/clusterbasis.c \
	Library/block.c \
	Library/h2matrix.c \
	Library/uniform.c \
	Library/kernelmatrix.c \
  	Library/coulomb.c \
	Library/lj.c \
	Library/tip4p.c \
	Library/rigid.c
	
SOURCES_libkips := \
	$(KIPS_CORE0) \
	$(KIPS_CORE1) \
	$(KIPS_CORE2)

EXTRA_HEADERS =

HEADERS_libkips := $(SOURCES_libkips:.c=.h) $(EXTRA_HEADERS)

OBJECTS_libkips := $(SOURCES_libkips:.c=.o)

DEPENDENCIES_libkips := $(SOURCES_libkips:.c=.d)

# ------------------------------------------------------------
# Test programs
# ------------------------------------------------------------

SOURCES_stable := \
	Tests/test_amatrix.c \
	Tests/test_interpolation.c \
	Tests/test_kernelmatrix.c \
  	Tests/test_coulomb.c \
	Tests/test_tip4p.c \
	Tests/test_eigenvalue.c

SOURCES_tests = \
	$(SOURCES_stable)

OBJECTS_tests := \
	$(SOURCES_tests:.c=.o)

DEPENDENCIES_tests := \
	$(SOURCES_tests:.c=.d)

PROGRAMS_tests := \
	$(SOURCES_tests:.c=)

# ------------------------------------------------------------
# Example programs
# ------------------------------------------------------------

SOURCES_examples =

OBJECTS_examples = \
	$(SOURCES_examples:.c=.o)

DEPENDENCIES_examples = \
	$(SOURCES_examples:.c=.d)

PROGRAMS_examples = \
	$(SOURCES_examples:.c=)

# ------------------------------------------------------------
# All files
# ------------------------------------------------------------

SOURCES := \
	$(SOURCES_libkips) \
	$(SOURCES_tests) \
	$(SOURCES_examples)

HEADERS := \
	$(HEADER_libkips)

OBJECTS := \
	$(OBJECTS_libkips) \
	$(OBJECTS_tests) \
	$(OBJECTS_examples)

DEPENDENCIES := \
	$(DEPENDENCIES_libkips) \
	$(DEPENDENCIES_tests) \
	$(DEPENDENCIES_examples)

PROGRAMS := \
	$(PROGRAMS_tests) \
	$(PROGRAMS_examples)

# ------------------------------------------------------------
# Standard target
# ------------------------------------------------------------

all: $(PROGRAMS_tests) $(PROGRAMS_examples)

# ------------------------------------------------------------
# Build configuration
# ------------------------------------------------------------

ifeq ($(wildcard options.inc),)
$(OBJECTS): options.inc.default
include options.inc.default
else
$(OBJECTS): options.inc
include options.inc
endif

# ------------------------------------------------------------
# System-dependent parameters (e.g., name of compiler)
# ------------------------------------------------------------

ifeq ($(wildcard system.inc),)
$(OBJECTS): system.inc.linux
include system.inc.linux
else
$(OBJECTS): system.inc
include system.inc
endif

# ------------------------------------------------------------
# Rules for test programs
# ------------------------------------------------------------

$(PROGRAMS_tests): %: %.o
ifdef BRIEF_OUTPUT
	@echo Linking $@
	@$(CC) $(LDFLAGS) $< -o $@ -lkips $(LIBS) 
else
	$(CC) $(LDFLAGS) $< -o $@ -lkips $(LIBS) 
endif

$(PROGRAMS_tests) $(PROGRAMS_tools): libkips.a

$(OBJECTS_tests): %.o: %.c
ifdef BRIEF_OUTPUT
	@echo Compiling $<
	@$(GCC) -MT $@ -MM -I Library $< > $(<:%.c=%.d)
	@$(CC) $(CFLAGS) -I Library -c $< -o $@
else
	@$(GCC) -MT $@ -MM -I Library $< > $(<:%.c=%.d)
	$(CC) $(CFLAGS) -I Library -c $< -o $@
endif

-include $(DEPENDENCIES_tests) $(DEPENDENCIES_tools)
$(OBJECTS_tests): Makefile

# ------------------------------------------------------------
# Rules for example programs
# ------------------------------------------------------------

$(PROGRAMS_examples): %: %.o
ifdef BRIEF_OUTPUT
	@echo Linking $@
	@$(CC) $(LDFLAGS) $< -o $@ -lkips $(LIBS) 
else
	$(CC) $(LDFLAGS) $< -o $@ -lkips $(LIBS) 
endif

$(PROGRAMS_examples): libkips.a

$(OBJECTS_examples): %.o: %.c
ifdef BRIEF_OUTPUT
	@echo Compiling $<
	@$(GCC) -MT $@ -MM -I Library $< > $(<:%.c=%.d)
	@$(CC) $(CFLAGS) -I Library -c $< -o $@
else
	@$(GCC) -MT $@ -MM -I Library $< > $(<:%.c=%.d)
	$(CC) $(CFLAGS) -I Library -c $< -o $@
endif

-include $(DEPENDENCIES_examples)
$(OBJECTS_examples): Makefile

# ------------------------------------------------------------
# Rules for the Doxygen documentation
# ------------------------------------------------------------

doc:
	doxygen Doc/Doxyfile

# ------------------------------------------------------------
# Rules for the main library
# ------------------------------------------------------------

libkips.a: $(OBJECTS_libkips)
ifdef BRIEF_OUTPUT
	@echo Building $@
	@$(AR) $(ARFLAGS) $@ $(OBJECTS_libkips)
else
	$(AR) $(ARFLAGS) $@ $(OBJECTS_libkips)
endif

$(OBJECTS_libkips): %.o: %.c
ifdef BRIEF_OUTPUT
	@echo Compiling $<
	@$(GCC) -MT $@ -MM $< > $(<:%.c=%.d)
	@$(CC) $(CFLAGS) -c $< -o $@
else
	@$(GCC) -MT $@ -MM $< > $(<:%.c=%.d)
	$(CC) $(CFLAGS) -c $< -o $@
endif

-include $(DEPENDENCIES_libkips)
$(OBJECTS_libkips): Makefile

# ------------------------------------------------------------
# Useful additions
# ------------------------------------------------------------

.PHONY: clean cleandoc programs indent

clean:
	$(RM) -f $(OBJECTS) $(DEPENDENCIES) $(PROGRAMS) libkips.a

cleandoc:
	$(RM) -rf Doc/html Doc/latex

indent:
	indent -bap -br -nce -cdw -npcs \
	  -di10 -nbc -brs -blf -i2 -lp \
	  -T amatrix -T pamatrix -T pcamatrix \
	  -T avector -T pavector -T pcavector \
	  -T cluster -T pcluster -T pccluster \
	  -T block -T pblock -T pcblock \
	  -T uniform -T puniform -T pcuniform \
	  -T h2matrix -T ph2matrix -T pch2matrix \
	  $(SOURCES)

coverage:
	mkdir Coverage > /dev/null 2>&1; \
	lcov --base-directory . --directory . --capture \
	--output-file Coverage/coverage.info && \
	genhtml -o Coverage Coverage/coverage.info

cleangcov:
	$(RM) -rf Library/*.gcov Library/*.gcda Library/*.gcno \
	Tests/*.gcov Tests/*.gcda Tests/*.gcno;
	$(RM) -rf Coverage
