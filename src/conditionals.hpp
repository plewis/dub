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
