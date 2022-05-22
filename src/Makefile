# SWM Makefile

include macros.make

ICOSIO_FC   = $(if $(filter $(PAR),serial),$(FCserial),$(FC))
ICOSIO_ARGS = CPP=$(CPP) CPPFLAGS="$(CPP_FLAGS) $(DEFINES) $(BUILD_DEFINES) -DNOGRIB" FC="$(ICOSIO_FC)" FFLAGS="$(ICOSIOFLAGS)"
SRCDIR      = $(CURDIR)

.PHONY: all build clean info

all: info build

info:
	$(info Inside top level Makefile)
	$(info ARCH   = $(ARCH))
	$(info HW     = $(HW))
	$(info PAR    = $(PAR))
	$(info MAXNZ  = $(MAXNZ))
	$(info SRCDIR = $(SRCDIR))
	$(info DOUBLE_PRECISION = $(DOUBLE_PRECISION))

build:
	mkdir -p bin lib     # YGYU 15-Apr-2020 bc -o bin/exe no longer work
	$(MAKE) $(GMAKEMINUSJ) $(ICOSIO_ARGS)               -C icosio
	$(MAKE) $(GMAKEMINUSJ) SRCDIR=$(SRCDIR)             -C utils
	$(MAKE) $(GMAKEMINUSJ)                              -C cntl
	$(MAKE) $(GMAKEMINUSJ)                              -C libsrc
	$(MAKE) $(GMAKEMINUSJ)                              -C dummygptl
	$(MAKE) SRCDIR=$(SRCDIR)                            -C swm
	$(MAKE)  					    -C stswm/src

clean:
	$(MAKE) -C icosio    clean
	$(MAKE) -C utils     clean
	$(MAKE) -C libsrc    clean
	$(MAKE) -C cntl      clean
	$(MAKE) -C dummygptl clean
	$(MAKE) -C swm       clean
	$(RM) bin/* lib/*
