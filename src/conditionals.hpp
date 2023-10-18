// Currently this is the only model implemented - do not uncomment!
#define USE_JUKE_CANTOR_MODEL
//#define DEBUGGING
//#define DEBUG_COAL_LIKE
//#define DEBUG_PARTIAL_STORE

// Use for saving continous parameters to a file for input into LoRaD software
// Note: this has very limited value. Can use it to check SMC marginal likelihoods
// only for the case where every gene tree has the same topology (and relative
// coalescence times!) and, even then, one can only compute p(D | S, theta, lambda),
// not p(D)
//#define SAVE_PARAMS_FOR_LORAD

// If defined, species are represented as bits set in an unsigned long
// If not defined, species are std::set
// Leaf nodes in gene trees are always assigned to one species, so
// their set would have one element or the unsigned long would have one bit set
// Ancestral species are unions of their descendant species.
#define SPECIES_IS_BITSET

//#define LOG_MEMORY
//#define USING_SIGNPOSTS
//#define USING_MPI
//#define ENABLE_PAUSE
#define PLOT_MULTIPLE_TRY_UPDATES

#if defined(USING_MPI)
#   include <mpi.h>
#endif

#if defined(ENABLE_PAUSE)
#   define DEBUG_PAUSE(s) cout << "paused: " << s << endl; cin >> dummy_char; cout << endl;
#else
#   define DEBUG_PAUSE(s)
#endif
