
EXECUTABLES = DPRM3.exe

DPRM3_OBJECTS   = angleMod.o ArrMod.o bdryMod.o   DPRM3Mod.o sspMod.o ReflectMod.o                       WriteRay.o              influence.o                           ReadEnvironmentBell.o DPRM3.o ../misc/libmisc.a
DPRM3_LIBS = -L../misc -lmisc

all: $(EXECUTABLES)
	@echo "DPRM3 built"
	@echo "*************"
	@echo " "

install: all
	for f in $(EXECUTABLES) ; do \
		echo "----- Installing $$f"; cp -p $$f ../bin; \
	done

clean:
	-rm -f *.o *.mod *.exe *_genmod.f90

# Declare implicit rules (for make utilities with missing or broken build rules)

.SUFFIXES: .a .o .mod .f .f90

.f.o:
	$(FC) -c $(FFLAGS) $< 

.f90.o:
	$(FC) -c $(FFLAGS) $< 

.f90.mod:
	$(FC) -c $(FFLAGS) $< 

# Declare explicit dependencies (and explicit build rules as needed)

DPRM3.exe: $(DPRM3_OBJECTS)
	$(FC) -o $@ $(FFLAGS) $(DPRM3_OBJECTS) $(DPRM3_LIBS)


angleMod.o anglemod.mod: angleMod.f90
	$(FC) -c $(FFLAGS) $< 

ArrMod.o arrmod.mod: ArrMod.f90
	$(FC) -c $(FFLAGS) $< 

bdryMod.o bdrymod.mod: bdryMod.f90
	$(FC) -c $(FFLAGS) $< 

bdry3DMod.o bdry3dmod.mod: bdry3DMod.f90
	$(FC) -c $(FFLAGS) $< 

DPRM3Mod.o DPRM3mod.mod: DPRM3Mod.f90
	$(FC) -c $(FFLAGS) $< 

Cone.o cone.mod: Cone.f90
	$(FC) -c $(FFLAGS) $<

influence.o influence.mod: influence.f90
	$(FC) -c $(FFLAGS) $< 

influence3D.o influence3d.mod: influence3D.f90
	$(FC) -c $(FFLAGS) $< 

RayNormals.o raynormals.mod: RayNormals.f90
	$(FC) -c $(FFLAGS) $<

ReadEnvironmentBell.o readenvironmentbell.mod: ReadEnvironmentBell.f90
	$(FC) -c $(FFLAGS) $< 

ReflectMod.o   reflectmod.mod:   ReflectMod.f90   Cone.o
	$(FC) -c $(FFLAGS) $<

Reflect3DMod.o reflect3dmod.mod: Reflect3DMod.f90 Cone.o
	$(FC) -c $(FFLAGS) $<

sspMod.o sspmod.mod: sspMod.f90
	$(FC) -c $(FFLAGS) $< 

Step3DMod.o step3dmod.mod: Step3DMod.f90
	$(FC) -c $(FFLAGS) $<

WriteRay.o writeray.mod: WriteRay.f90
	$(FC) -c $(FFLAGS) $<

DPRM3.o:   arrmod.mod anglemod.mod bdrymod.mod   readenvironmentbell.mod DPRM3mod.mod sspmod.mod reflectmod.mod                  writeray.mod raynormals.mod influence.mod                 DPRM3.f90


