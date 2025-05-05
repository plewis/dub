// Currently this is the only model implemented - do not uncomment!
#define USE_JUKE_CANTOR_MODEL

// Alters output to be more helpful to developers when debugging
// (but unhelpful to users, however)
//#define OUTPUT_FOR_DEBUGGING

// Calculates log likelihood and subtracts previous log likelihood to
// check whether log weight calculation is correct.
// Should be undefined in production version for speed.
#define DEBUG_CHECK_WEIGHTS

//#define ADHOC_REMOVE_BEFORE_RELEASE
//#define DEBUG_COALLIKE

//#define DEBUGGING_INITFROMPARTICLE
//#define HACK_FOR_SNAKE_ATP

//#define RANDOM_LOCUS_ORDERING

//#define FOSSILS

#define REUSE_PARTIALS
//#define USE_HEATING
//#define SYSTEMATIC_FILTERING

#define LAZY_COPYING
#define USING_MULTITHREADING

#define OUTPUT_DEBUG(a)

