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

#if defined(WITH_OPENMP)
#define __TotalNrOfThreads      OMP_GET_NUM_THREADS()
#define __CurrentThreadNum      OMP_GET_THREAD_NUM() + 1
#else
#define __TotalNrOfThreads      1
#define __CurrentThreadNum      1
#endif

   /* These are modules that are used in most parts of the code */
   USE ErrorTrap
   USE MyConsts
#if defined(WITH_OPENMP)
   USE omp_lib
#endif
   