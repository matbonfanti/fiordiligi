/* This is the general preprocess include file. It should be included in every fortran code at
   the place where the USE statements are located. This ensures the correct working of all the 
   macros and defines located here and supplied by the make files... 
*/

/*==============================================================================*
 * Define alias for log file management.                                        *
 *==============================================================================*/

#if defined(LOG_FILE) 
#define __LOG_UNIT              19
#define __OPEN_LOG_FILE         OPEN( UNIT=19, FILE=LOG_FILE, ACTION="write", POSITION="append" )
#define __CLOSE_LOG_FILE        CLOSE( UNIT=19 )
#endif

/*==============================================================================*
 * Define/undefine printing of graphene slab for vhsticking 3D function.        *
 *==============================================================================*/

#undef __PRINT_VHSTICKING_FUNCTION
#if defined(__PRINT_VHSTICKING_FUNCTION) 
#define __VHSTICKING_UNIT              55
#define __OPEN_VHSTICKING_FILE         OPEN( UNIT=55, FILE="vhsticking.f90", ACTION="write" )
#define __CLOSE_VHSTICKING_FILE        CLOSE( UNIT=55 )
#endif

/*==============================================================================*
 * Define/undefine printing of spectral density.                                *
 *==============================================================================*/

#undef __PRINT_SPECTRAL_DENSITY

/*==============================================================================*
 * Define/undefine ZH dependent coupling for scattering simulations.            *
 *==============================================================================*/

#undef __ZH_DEPENDENT_COUPLING

/*==============================================================================*
 * Define alias for OMP parallelization.                                        *
 *==============================================================================*/

#if defined(WITH_OPENMP)
#define __OMP_TotalNrOfThreads      OMP_GET_NUM_THREADS()
#define __OMP_CurrentThreadNum      OMP_GET_THREAD_NUM() + 1
#define __OMP_OnlyMasterBEGIN       IF ( OMP_GET_THREAD_NUM() == 0 ) THEN;
#define __OMP_OnlyMasterEND         ; END IF
#else
#define __OMP_TotalNrOfThreads      1
#define __OMP_CurrentThreadNum      1
#define __OMP_OnlyMasterBEGIN   
#define __OMP_OnlyMasterEND     
#endif

/*==============================================================================*
 * Load common modules.                                                         *
 *==============================================================================*/

   /* These are modules that are used in most parts of the code */
   USE ErrorTrap
   USE MyConsts
   USE MyLinearAlgebra
#if defined(WITH_OPENMP)
   USE omp_lib
#endif
   
