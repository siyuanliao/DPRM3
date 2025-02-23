
AR = ar

# Intel compiler needed the following definition for the archiver:
# AR = xiar

ARFLAGS = r

RANLIB = ranlib

# include norm2.o below, if your compiler is missing that intrinsic
OBJECTS = ReadHeader.o beampattern.o MathConstants.o RefCoef.o SourceReceiverPositions.o  pchipMod.o AttenMod.o sspMod.o WriteHeader.o errout.o interpolation.o MergeVectors.o munk.o ReadEnvironmentMod.o SortMod.o splinec.o subtabulate.o twersky.o calculateweights.o xerror.o norms.o cross_products.o monotonicMod.o PolyMod.o PekerisRoot.o

all: libmisc.a

install: all

clean:
	-rm -f *.a *.mod *.o *_genmod.f90

# Declare implicit rules (for make utilities with missing or broken build rules)

.SUFFIXES: .a .o .mod .f .f90

.f.o:
	$(FC) -c $(FFLAGS) $< 

.f90.o:
	$(FC) -c $(FFLAGS) $< 

.f90.mod:
	$(FC) -c $(FFLAGS) $< 


# Declare explicit dependencies (and explicit build rules as needed)

libmisc.a: $(OBJECTS)
	$(AR) $(ARFLAGS) $@ $(OBJECTS)
	$(RANLIB) $@
	@echo "Routines needed by SCOOTER\SPARC\BELLHOP\BOUNCE built"

subtabultate.o subtabulate.mod: subtabulate.f90
	$(FC) -c $(FFLAGS) $<

calculateweights.o calculateweights.mod: calculateweights.f90
	$(FC) -c $(FFLAGS) $<

SourceReceiverPositions.o sourcereceiverpositions.mod: SourceReceiverPositions.f90 monotonicmod.mod sortmod.mod subtabulate.mod
	$(FC) -c $(FFLAGS) $< 

ReadHeader.o: sourcereceiverpositions.mod ReadHeader.f90

PolyMod.o polymod.mod: PolyMod.f90
	$(FC) -c $(FFLAGS) $< 

AttenMod.o: AttenMod.f90 MathConstants.mod
	$(FC) -c $(FFLAGS) $<

beampattern.o beampattern.mod: beampattern.f90
	$(FC) -c $(FFLAGS) $<

interpolation.o interpolation.mod: interpolation.f90
	$(FC) -c $(FFLAGS) $<

MathConstants.o mathconstants.mod: MathConstants.f90
	$(FC) -c $(FFLAGS) $< 

RefCoef.o refcoef.mod: RefCoef.f90 polymod.mod
	$(FC) -c $(FFLAGS) $< 

sspMod.o sspmod.mod: sspMod.f90 pchipmod.mod attenmod.mod
	$(FC) -c $(FFLAGS) $< 

pchipMod.o pchipmod.mod: pchipMod.f90
	$(FC) -c $(FFLAGS) $< 

monotonicMod.o monotonicmod.mod: monotonicMod.f90
	$(FC) -c $(FFLAGS) $< 

SortMod.o sortmod.mod: SortMod.f90
	$(FC) -c $(FFLAGS) $< 

PekerisRoot.o: PekerisRoot.f90

WriteHeader.o: WriteHeader.f90 sourcereceiverpositions.mod
errout.o: errout.f90
MergeVectors.o: MergeVectors.f90
munk.o: munk.f90
norms.o norms.mod: norms.f90
	$(FC) -c $(FFLAGS) $< 


ReadEnvironmentMod.o readenvironmentmod.mod: ReadEnvironmentMod.f90 sspmod.mod
splinec.o: splinec.f90
twersky.o: twersky.f90
xerror.o: xerror.f90
cross_products.o cross_products.mod: cross_products.f90
	$(FC) -c $(FFLAGS) $< 

