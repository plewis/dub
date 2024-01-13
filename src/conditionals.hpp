// Currently this is the only model implemented - do not uncomment!
#define USE_JUKE_CANTOR_MODEL

//temporary!
//#define POLTMP1

// Uncomment for release version
//#define DEBUGGING_SANITY_CHECK

// Uncomment to produce memory report for partials and data
//#define LOG_MEMORY

// Uncomment to use signposts in Instruments
//#define USING_SIGNPOSTS

// Uncomment to use Teh's prior-post method for coalescent events
//#define PRIOR_POST

//#define USING_MPI
//#define USING_MULTITHREADING

//#define MINIMIZE_PARTIALS

#if defined(USING_MPI)
#   error MPI version not yet completed
#   include <mpi.h>
#endif

#if defined(USING_MPI) && defined(USING_MULTITHREADING)
#error USING_MPI and USING_MULTITHREADING cannot both be defined at the same time
#endif
