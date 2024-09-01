//#define NDEBUG
#include <cassert>
#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <memory>   // shared_ptr
#include <new>      // used by Mallocator: bad_alloc, bad_array_new_length
#include <numeric>
#include <queue>
#include <regex>
#include <set>
#include <stack>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/program_options.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include "ncl/nxsmultiformat.h"

using namespace std;
using boost::format;

// See https://www.jviotti.com/2022/02/21/emitting-signposts-to-instruments-on-macos-using-cpp.html
#if defined(USING_SIGNPOSTS)
#   include <os/log.h>
#   include <os/signpost.h>
    os_log_t log_handle;
    os_signpost_id_t signpost_id;
#endif

#include "conditionals.hpp"
#include "xproj.hpp"
#include "lot.hpp"
#include "smcglobal.hpp"
#include "stopwatch.hpp"
#include "genetic-code.hpp"
#include "datatype.hpp"
#include "partition.hpp"
#include "data.hpp"
#include "split.hpp"
#include "partial_store.hpp"
#include "node.hpp"
#include "forest.hpp"
#include "species-forest.hpp"
#include "policy-parallel-none.hpp"
#include "smc.hpp"
#include "particle.hpp"
#include "gene-forest.hpp"
#include "particle-func.hpp"
#include "smc-func.hpp"
//#include "policy-parallel-mt.hpp"
//#include "policy-parallel-mpi.hpp"

#if defined(DLIB_EXPERIMENT)
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
#endif

#include "proj.hpp"

using namespace proj;

#if defined(USING_MPI) && defined(USING_MULTITHREADING)
#error USING_MPI and USING_MULTITHREADING are mutually exclusive: undefine one or the other to compile
#endif

#if defined(USING_MPI)
    int my_rank = 0;
    int ntasks = 0;

    unsigned my_first_particle = 0;
    unsigned my_last_particle = 0;

    void output(format & fmt, unsigned level) {
        if (my_rank == 0 && G::_verbosity > 0 && level <= G::_verbosity)
            cout << str(fmt);
    }

    void output(string msg, unsigned level) {
        if (my_rank == 0 && G::_verbosity > 0 && level <= G::_verbosity)
            cout << msg;
    }
#else
    void output(format & fmt, unsigned level) {
        if (G::_verbosity > 0 && level <= G::_verbosity)
            cout << str(fmt);
    }

    void output(string msg, unsigned level) {
        if (G::_verbosity > 0 && level <= G::_verbosity)
            cout << msg;
    }
#endif

#if defined(LOG_MEMORY)
    char dummy_char;
    ofstream memfile;
#endif

//POLWAS Lot                                 rng;
Lot::SharedPtr                      rng(new Lot());
StopWatch                           stopwatch;
PartialStore                        ps;

#if defined(LOG_MEMORY)
vector<unsigned>                    Partial::_nconstructed;
vector<unsigned>                    Partial::_ndestroyed;
vector<unsigned>                    Partial::_max_in_use;
vector<unsigned long>               Partial::_bytes_per_partial;
unsigned long                       Partial::_total_max_in_use  = 0;
unsigned long                       Partial::_total_max_bytes   = 0;
unsigned                            Partial::_nstates           = 4;
#endif

string                              G::_species_tree_ref_file_name = "";
string                              G::_gene_trees_ref_file_name = "";

string                              G::_debugging_text          = "doof.txt";
bool                                G::_debugging               = false;

unsigned                            G::_nthreads        = 1;
#if defined(USING_MULTITHREADING)
mutex                               G::_mutex;
mutex                               G::_gene_forest_clear_mutex;
mutex                               G::_debug_mutex;
//vector<unsigned>                    G::_thread_first_gene;
//vector<unsigned>                    G::_thread_last_gene;
#endif

// Let Nup be the number of unique particles
// Let Np be the nominal number of particles (i.e. count = 1 for each particle)
// Nup < Np if count > 1 for some particles
// 0: Save all species trees even if they are identical to others
//    File will contain Np species trees
// 1: Save species tree from each unique particle and, for each, show that particle's count
//    File will contain Nup species trees
// 2: Save unique species trees along with their frequency in the sample
//    Is this effectively the same as compression level 1?
//    Creates map in which species tree newicks are keys and cumulative count is value,
//    but will there ever be exactly the same newick in different particles?
unsigned                            G::_treefile_compression = 0;

unsigned                            G::_verbosity          = 3;

unsigned                            G::_nstates            = 4;

unsigned                            G::_ntaxa              = 0;
vector<string>                      G::_taxon_names;
map<string, unsigned>               G::_taxon_to_species;

unsigned                            G::_nspecies           = 0;
G::species_t                        G::_species_mask       = (G::species_t)0;
vector<string>                      G::_species_names;
map<unsigned,unsigned>              G::_nexus_taxon_map;

unsigned                            G::_ngenes              = 0;
vector<string>                      G::_gene_names;
vector<unsigned>                    G::_nsites_per_gene;
map<unsigned, double>               G::_relrate_for_gene;

double                              G::_phi                 = 1.0;
double                              G::_theta               = 0.05;
double                              G::_lambda              = 1.0;

double                              G::_invgamma_shape      = 2.0;
bool                                G::_theta_mean_frozen   = false;
double                              G::_theta_mean_fixed    = -1.0;
double                              G::_theta_prior_mean    = 1.0;
double                              G::_theta_proposal_mean = 0.1;
double                              G::_lambda_prior_mean   = 1.0;

//bool                                G::_update_theta       = true;
//bool                                G::_update_lambda      = true;

double                              G::_small_enough         = 0.00001;

unsigned                            G::_nparticles           = 500;
unsigned                            G::_nkept                = 500;
unsigned                            G::_nparticles2          = 1000;

//temporary!
map<unsigned, vector<G::SpecLog> >  G::_speclog;

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required in order to use infinity()");
double                              G::_infinity = numeric_limits<double>::infinity();
double                              G::_negative_infinity = -numeric_limits<double>::infinity();

bool                                G::_prior_post        = false;

PartialStore::leaf_partials_t       GeneForest::_leaf_partials;

const double                        Node::_smallest_edge_length = 1.0e-12;

string                              Proj::_program_name         = "smc";
unsigned                            Proj::_major_version        = 3;
unsigned                            Proj::_minor_version        = 0;

GeneticCode::genetic_code_definitions_t GeneticCode::_definitions = { 
                             // codon order is alphabetical: i.e. AAA, AAC, AAG, AAT, ACA, ..., TTT
    {"standard",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"vertmito",             "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"yeastmito",            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"moldmito",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"invertmito",           "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"ciliate",              "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF"},
    {"echinomito",           "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"euplotid",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF"},
    {"plantplastid",         "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"altyeast",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"ascidianmito",         "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"altflatwormmito",      "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF"},
    {"blepharismamacro",     "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF"},
    {"chlorophyceanmito",    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF"},
    {"trematodemito",        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"scenedesmusmito",      "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF"},
    {"thraustochytriummito", "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF"}
};

int main(int argc, const char * argv[]) {

#if defined(USING_MPI)
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
#   if defined(LOG_MEMORY)
        memfile.open(str(format("allocs-%d.txt") % my_rank));
#   endif
#elif defined(LOG_MEMORY)
    memfile.open("allocs.txt");
#endif

#if defined(USING_SIGNPOSTS)
    log_handle  = os_log_create("edu.uconn.eeb.phylogeny", OS_LOG_CATEGORY_POINTS_OF_INTEREST);
    signpost_id = os_signpost_id_generate(log_handle);
    assert(signpost_id != OS_SIGNPOST_ID_INVALID);
#endif

    Proj proj;
    bool normal_termination = true;
    try {
        proj.processCommandLineOptions(argc, argv);
        
        StopWatch sw;
        sw.start();
        proj.run();
        double total_seconds = sw.stop();
        output(format("\nTotal time: %.5f seconds\n") % total_seconds, 1);
    }
    catch(std::exception & x) {
        cerr << str(format("Exception: %s\n") % x.what());
        cerr << "Aborted.\n";
        normal_termination = false;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        normal_termination = false;
    }
    
#if defined(LOG_MEMORY)
    if (normal_termination) {
        proj.memoryReport(memfile);
        ps.memoryReport(memfile);
    }
    memfile.close();
#endif

    return 0;
}
