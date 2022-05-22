# dynamics Makefile (for NIM dynamics compiled for the CPU)

include ../macros.make
include NIM_DYNAMICS_OBJS

.SMS     = $(SMS)/bin/ppp modinfo
CNTLOBJS = $(addsuffix .o, $(basename $(CNTLSRCS)))
CNTLSRCS = $(wildcard $(CNTLDIR)/*.F90)
BASICOBJS= ../libsrc/basiclib/mathut.o   ../lib/libslint.a  ../lib/libmath.a
LIBS     = $(GPTLLIB) $(INCNETCDF) $(LIBNETCDF) -L../lib -lmath -lslint

.PHONY: clean objects1 objects2

%.o: %.f
	$(FCX) -c $<

all: DEPENDENCIES objects1 objects2 $(SWMEXE)

objects1:
	$(MAKE) $(GMAKEMINUSJ) -f Makefile.serial.cpu SFLAGS="$(SFLAGS)" $(SEROBJS)
ifeq ($(HW),gpu)
	/bin/cp cuda/*.cu .
	$(MAKE) $(GMAKEMINUSJ) -f Makefile.serial.gpu $(COMOBJS) $(GPUOBJS)
  LIBS += $(LIBCUDA)
else
	$(MAKE) $(GMAKEMINUSJ) -f Makefile.serial.cpu SFLAGS="" $(COMOBJS) $(CPUOBJS)
endif

objects2:
ifeq ($(PAR),serial)
	$(MAKE) $(GMAKEMINUSJ) -f Makefile.serial.cpu SFLAGS="" $(SMSOBJS)
else
	$(.SMS) ../physics/grims/common/rdparm8.f90   . # .sms file for grims_setup
	$(.SMS) ../physics/grims/common/comconsts.f90 . # .sms file for grims_setup
	$(MAKE) $(GMAKEMINUSJ) -f Makefile.sms $(SMSOBJS)
  LIBS += -L$(SMS)/lib -l$(SMSLIBNAME)
endif

ifeq ($(HW),gpu)
  OBJS = $(CNTLOBJS) $(SMSOBJS) $(SEROBJS) $(COMOBJS) $(GPUOBJS)
else
  OBJS = $(CNTLOBJS) $(SMSOBJS) $(SEROBJS) $(COMOBJS) $(CPUOBJS)
endif


$(SWMEXE): $(OBJS) $(ICOSIOOBJS) $(UTILSOBJS) $(GPTLOBJS) $(BASICOBJS)
	$(FCX) -o $@ $(MODDIRS) $(OMPFLAG) $(MICFLAGS) $(OBJS) \
        $(ICOSIOOBJS) $(UTILSOBJS) $(GPTLOBJS) $(LIBS)


DEPENDENCIES:
	$(RM) Filepath Srcfiles
	echo "." > Filepath
#	ls *.F90 *.f > Srcfiles
	ls *.F90  > Srcfiles
	$(MKDEPENDS) -m Filepath Srcfiles > $@

clean:
	$(RM) DEPENDENCIES .ACC* *.cu *.f90 *.m4 *.mod *.o *.sms
