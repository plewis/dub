// Currently this is the only model implemented - do not uncomment!
#define USE_JUKE_CANTOR_MODEL

// //temporary!
//#define POLTMP

// Comment out for release version
//#define DEBUGGING_SANITY_CHECK

// Comment out for release version
//#define DEBUG_SECOND_LEVEL

// Uncomment to let each primary particle have its own mean theta
// Secondary particles use their inherited mean theta to
// integrate out individual species thetas using Jones (2017)
#define EST_THETA

// Uncomment to produce memory report for partials and data
//#define LOG_MEMORY

// Uncomment to use signposts in Instruments
//#define USING_SIGNPOSTS

// Uncomment to use Teh's prior-post method for coalescent events
//#define PRIOR_POST

//#define MINIMIZE_PARTIALS

#if defined(USING_MPI)
#   error MPI version not yet completed
#   include <mpi.h>
#endif
