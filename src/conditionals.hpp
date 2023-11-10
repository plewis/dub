// Currently this is the only model implemented - do not uncomment!
#define USE_JUKE_CANTOR_MODEL

// Uncomment to produce memory report for partials and data
#define LOG_MEMORY

// Uncomment to sacrifice speed for lower memory requirements
//TODO: actually much slower and uses much more memory if LOW_MEM is defined!
//  need to limit how many partials are stored for each gene?
//  On test/simulated allocated 4362935 partials under LOW_MEM vs 88968 otherwise
//  partials allocated lomem/himem: 4362935/88968 = 49.0
//  real: 24.711/9.489 = 2.6
//  user:18.585/8.574  = 2.2
//  sys:6.039/0.827    = 7.3
//#define LOW_MEM

// Uncomment to use signposts in Instruments
//#define USING_SIGNPOSTS

// Uncomment to use Teh's prior-post method for coalescent events
//#define PRIOR_POST

//#define USING_MPI

#define USING_MULTITHREADING
