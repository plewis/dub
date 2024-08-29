// Currently this is the only model implemented - do not uncomment!
#define USE_JUKE_CANTOR_MODEL

//#define POLTMP

// Comment out for release version
//#define DEBUGGING_SANITY_CHECK

// Comment out for release version
//#define DEBUG_SECOND_LEVEL

// Comment out for release version
// If defined, recalculates coalescent likelihood before proposing
// new increment/join to be sure _prev_log_coallike is correct
//#define CHECK_SECOND_LEVEL_PREV_LOG_LIKELIHOOD

// Uncomment to sanity check log weight for particle after Particle::proposeCoalescence
//#define DEBUG_CHECK_LOGWEIGHT

// Uncomment to let each primary particle have its own mean theta
// Secondary particles use their inherited mean theta to
// integrate out individual species thetas using Jones (2017)
#define EST_THETA

// Experimental
//#define WEIGHT_BY_GENE_LENGTH

// Experimental: devote each step to advancing a single locus
#define ONE_LOCUS_PER_STEP

// Uncomment to use simpler specification of theta fixedness
//#define SIMPLIFY_THETA

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

#if defined(MINIMIZE_PARTIALS) && defined(ONE_LOCUS_PER_STEP)
#   error cannot currently combine MINIMIZE_PARTIALS and ONE_LOCUS_PER_STEP
#endif
