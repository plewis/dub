// Currently this is the only model implemented - do not uncomment!
#define USE_JUKE_CANTOR_MODEL

// Alters output to be more helpful to developers when debugging
// (but unhelpful to users, however)
//#define OUTPUT_FOR_DEBUGGING

// Calculate likelihood using the existing forest with the rest of
// the tree filled in by UPGMA.
#define UPGMA_WEIGHTS
//#define DEBUG_UPGMA

//#define ADHOC_REMOVE_BEFORE_RELEASE
//#define DEBUG_COALLIKE

//#define DEBUGGING_INITFROMPARTICLE
//#define HACK_FOR_SNAKE_ATP

//#define RANDOM_LOCUS_ORDERING

//#define FOSSILS

#define REUSE_PARTIALS
//#define USE_HEATING
#define PRECALC_JC_TRANSITION_PROBS
//#define SYSTEMATIC_FILTERING

//#define OUTPUT_DEBUG(a) output((a), G::LogCateg::INFO)
#define OUTPUT_DEBUG(a)
