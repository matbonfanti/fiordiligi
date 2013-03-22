
# ******************** MethaneOnSurfaceCRP MAKEFILE **************** 

#----------------------------------------------------------------------------
#                         USER DEFINABLE OPTIONS
#----------------------------------------------------------------------------

# Executable name
EXENAME = JK6_v2

# Print relevant output to log file
LOGFILE = yes
LOGNAME = jk6_v2.log

# Compiler ( gfortran, ifort )
FC = ifort

# Debugging options ( yes or no )
DEBUG = no  

# Optimization level
OPTLEVEL = 3

# linking LAPACK and BLAS 
LAPACK = no

# Intel compiler version 
# required only if FC=ifort and LAPACK=yes
INTELVERS = 11

# Compile with standard real 8 (see details about the flags for each compiler...)
REAL8 = yes

#----------------------------------------------------------------------------
#                             STRIP ALL SPACES
#----------------------------------------------------------------------------

# Strip leading and trailing spaces from all variables.
FC := $(strip ${FC})
DEBUG := $(strip ${DEBUG})
OPTLEVEL := $(strip ${OPTLEVEL})
LAPACK := $(strip ${LAPACK})
INTELVERS := $(strip ${INTELVERS})
LOGFILE := $(strip ${LOGFILE})
LOGNAME := $(strip ${LOGNAME})



#----------------------------------------------------------------------------
#                       Compiler specific statements.
#----------------------------------------------------------------------------

ifeq (${FC},gfortran)

   # Optimization flag
   O0FLAGS  = -O0
   O1FLAGS  = -O1
   O2FLAGS  = -O2
   O3FLAGS  = -O3

   # Debug flag
   DEBUGFLG = -g -fbounds-check

   # GNU Lapack and Blas flags 
   LAPACKFLG = -llapack -lblas

   # ATLAS Lapack and Blas flags
   LAPACKFLG =  -L/usr/lib/atlas/ -llapack -lf77blas -lcblas -latlas

   # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
	DATAFLG = -fdefault-real-8
   endif

   # Flag to specify the position of mod files
   MODULEFLG = -fintrinsic-modules-path

endif

ifeq (${FC},ifort)

   # Optimization flags
   O0FLAGS  = -O0
   O1FLAGS  = -O1
   O2FLAGS  = -O2
   O3FLAGS  = -O3

   # Debug flags
   DEBUGFLG  =  -g -traceback   -fpe-all=0    -debug all -check all 

   # MKL flags
   ifeq (${INTELVERS},11)
      LAPACKFLG = -lmkl_lapack  -lmkl_intel -lmkl_sequential -lmkl_core -lguide -lpthread
   endif
   ifeq (${INTELVERS},10)
      LAPACKFLG = -lmkl -lmkl_lapack -lguide -lpthread 
   endif

   # Data type
   DATAFLG =
   ifeq (${REAL8},yes)
	DATAFLG = -r8 -i4
   endif
 
   # Flag to specify the position of mod files
   MODULEFLG = -I

endif


#----------------------------------------------------------------------------
#              Setup preprocessing options 
#----------------------------------------------------------------------------

PPDEFINE = 

# Setup preprocessing debug options
ifeq (${DEBUG}, yes)
   PPDEFINE += -DVERBOSE_OUTPUT
endif

# Preprocess with lapack calls
ifeq (${LAPACK}, yes)
   PPDEFINE += -DWITH_LAPACK
endif


#----------------------------------------------------------------------------
#              Setup linking and compilation flags
#----------------------------------------------------------------------------

# initialize flags
COMPILEFLG = 
LINKFLG    = 
LIBFLG     = 

# if debugging set the appropriate flags
ifeq (${DEBUG}, yes)
   COMPILEFLG += ${DEBUGFLG}
endif

# Set flags for defining standard variable kinds
COMPILEFLG +=  ${DATAFLG} 

# If lapack, add the linking options
ifeq (${LAPACK}, yes)
   LIBFLG += ${LAPACKFLG}
endif


#----------------------------------------------------------------------------
#              Determine the optimization level to be used.
#----------------------------------------------------------------------------

# if debugging override input optimization level
ifeq (${DEBUG}, yes)
   OPTLEVEL = 0
endif

# Define optimization level
OPTFLAGS       = ${O0FLAGS}
ifeq (${OPTLEVEL},1)
  OPTFLAGS       = ${O1FLAGS}
endif
ifeq (${OPTLEVEL},2)
  OPTFLAGS       = ${O2FLAGS}
endif
ifeq (${OPTLEVEL},3)
  OPTFLAGS       = ${O3FLAGS}
endif

COMPILEFLG += ${OPTFLAGS}
#LINKFLG    += ${OPTFLAGS}


#----------------------------------------------------------------------------
#                      List of directories
#----------------------------------------------------------------------------

SRCDIR  = Source
PPDIR   = PreProcessed
OBJDIR  = Objects
TESTDIR = Tests/Source
EXEDIR  = Executables

#----------------------------------------------------------------------------
#                      List of object files
#----------------------------------------------------------------------------

# Define list of object from the list of all f90 files in the directory
OBJSWITHMAIN =$(patsubst Source/%,Objects/%,$(patsubst %.f90,%.o,$(wildcard ${SRCDIR}/*.f90)))
OBJS =$(patsubst %/Main.o,,${OBJSWITHMAIN})

#----------------------------------------------------------------------------
#       Construct the compile, link, preprocess variables.
#----------------------------------------------------------------------------

# Compile command: ${COMPILE} <source>
COMPILE                 = ${FC} ${COMPILEFLG} ${MODULEFLG} ${OBJDIR} -c 

# Link command: ${LINK} <exe name> <objects> <libflags>
LINK                    = ${FC} ${LINKFLG} ${MODULEFLG} ${OBJDIR} -o

# Preprocess: ${PREPROCESS} <to preprocess> <preprocessed>
PREPROCESS              = cpp -P -traditional ${PPDEFINE}

# Build static library : ${AR} <libraryname>.a <objects>
AR 			= ar cr

#----------------------------------------------------------------------------
#                       START OF MAKE RULES
#----------------------------------------------------------------------------

# Link objects to the produce the executable file 
${EXENAME} : ${SRCDIR}/Main.f90 ${OBJS} 
	${PREPROCESS} ${SRCDIR}/Main.f90 ${PPDIR}/Main.f90
	${COMPILE} ${PPDIR}/Main.f90
	${LINK} ${EXEDIR}/$@ Main.o $(OBJS) ${LIBFLG}
	rm Main.o

# Make target to build all the object files and assemble them
all : ${OBJS}

# Make a target object file by preprocessing and compiling the fortran code
${OBJDIR}/%.o : ${SRCDIR}/%.f90
	${PREPROCESS} ${SRCDIR}/$*.f90 ${PPDIR}/$*.f90
	${COMPILE} ${PPDIR}/$*.f90
	cp -p $*.o $(shell echo $* | tr A-Z a-z).mod ${OBJDIR}
	rm $*.o $(shell echo $* | tr A-Z a-z).mod

# Make target to build required directories
directories : ${PPDIR} ${OBJDIR} ${EXEDIR}
	mkdir -p ${PPDIR} ${OBJDIR} ${EXEDIR}

# Make documentation with doxygen
doc :
	doxygen Documentation/Doxyfile

# Remove compiled objects and related stuff
clean :
	rm -fr ${OBJDIR}/* ${PPDIR}/* 

# Clean documentation
clean-doc :
	rm -fr Documentation/html
	rm -fr Documentation/latex

# --------------------------------------------------------------------------------------------
# ---------------------------------- START WITH DEPENDENCIES NOW -----------------------------
# --------------------------------------------------------------------------------------------

# Very basic files, which everything depends on
COMMONDEP = Makefile ${OBJDIR}/ErrorTrap.o  ${OBJDIR}/MyConsts.o

# Set error and warning procedures
${OBJDIR}/ErrorTrap.o        : ${SRCDIR}/ErrorTrap.f90

# Define common physical and mathematical constants
${OBJDIR}/MyConsts.o         : ${SRCDIR}/MyConsts.f90 ${OBJDIR}/ErrorTrap.o

# Set procedure for reading quasi-free format input file
${OBJDIR}/InputField.o       : ${SRCDIR}/InputField.f90 ${COMMONDEP}

# Module containing the common data
${OBJDIR}/CommonData.o       : ${SRCDIR}/CommonData.f90 ${COMMONDEP}

# Module containing the potential energy surface
${OBJDIR}/PotentialModule.o  : ${SRCDIR}/PotentialModule.f90 ${OBJDIR}/CommonData.o ${COMMONDEP}

# Module containing the random number generator
${OBJDIR}/RandomNumberGenerator.o  : ${SRCDIR}/RandomNumberGenerator.f90 ${COMMONDEP}

# Module containing the integrator for the classical eq of motion
${OBJDIR}/ClassicalEqMotion.o  : ${SRCDIR}/ClassicalEqMotion.f90 ${OBJDIR}/RandomNumberGenerator.o ${COMMONDEP}

# Module containing the definitions of the independent oscillator model
${OBJDIR}/IndependentOscillatorsModel.o  : ${SRCDIR}/IndependentOscillatorsModel.f90 ${OBJDIR}/MyConsts.o ${OBJDIR}/PotentialModule.o \
                                           ${OBJDIR}/SplineInterpolator.o ${COMMONDEP}

# Module containing the spline interpolation subroutines
${OBJDIR}/SplineInterpolator.o : ${SRCDIR}/SplineInterpolator.f90 ${OBJDIR}/NRUtility.o ${COMMONDEP}

# Module containing the subroutine to write data in vtk format
${OBJDIR}/PrintTools.o : ${SRCDIR}/PrintTools.f90 ${COMMONDEP}

