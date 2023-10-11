//#define NDEBUG
#include <cassert>
#include <algorithm>
#include <climits>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
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

#include "ncl/nxsmultiformat.h"

using namespace std;
using boost::format;

#include "conditionals.hpp"

// See https://www.jviotti.com/2022/02/21/emitting-signposts-to-instruments-on-macos-using-cpp.html
#if defined(USING_SIGNPOSTS)
#   include <os/log.h>
#   include <os/signpost.h>
    os_log_t log_handle;
    os_signpost_id_t signpost_id;
#endif

#include "mallocator.hpp"
#include "stopwatch.hpp"
#include "xproj.hpp"
#include "lot.hpp"
#include "genetic-code.hpp"
#include "datatype.hpp"
#include "partition.hpp"
#include "data.hpp"
#include "split.hpp"
#include "partial_store.hpp"
#include "node.hpp"
#include "epoch.hpp"
#include "forest.hpp"
#include "gene-forest.hpp"
#include "species-forest.hpp"
#include "particle.hpp"
#include "epoch.hpp"
#include "proj.hpp"

using namespace proj;

unsigned Proj::_verbosity = 0;

#if defined(USING_MPI)
int my_rank = 0;
int ntasks = 0;

unsigned my_first_gene = 0;
unsigned my_last_gene = 0;

void output(format & fmt, unsigned level) {
    if (my_rank == 0 && Proj::_verbosity > 0 && level <= Proj::_verbosity) {
        cout << str(fmt);
    }
}

void output(string msg, unsigned level) {
    if (my_rank == 0 && Proj::_verbosity > 0 && level <= Proj::_verbosity) {
        cout << msg;
    }
}
#else
void output(format & fmt, unsigned level) {
    if (Proj::_verbosity > 0 && level <= Proj::_verbosity)
        cout << str(fmt);
}

void output(string msg, unsigned level) {
    if (Proj::_verbosity > 0 && level <= Proj::_verbosity)
        cout << msg;
}
#endif

#if defined(LOG_MEMORY)
    char dummy_char;
    ofstream memfile;
#endif

Lot                                 rng;
StopWatch                           stopwatch;
PartialStore                        ps;

#if defined(LOG_MEMORY)
vector<unsigned>                    Partial::_nconstructed;
vector<unsigned>                    Partial::_ndestroyed;
vector<unsigned>                    Partial::_max_in_use;
vector<unsigned long>               Partial::_bytes_per_partial;
unsigned                            Partial::_total_max_in_use  = 0;
unsigned                            Partial::_total_max_bytes   = 0;
unsigned                            Partial::_nstates           = 4;

unsigned                            Epoch::_nconstructed        = 0;
unsigned                            Epoch::_ndestroyed          = 0;
unsigned                            Epoch::_max_in_use          = 0;
unsigned                            Epoch::_bytes_base_epoch = (unsigned long)(
    2*sizeof(int)               // _type, _gene
    + sizeof(double)            // _height
    + sizeof(bool)              // _valid
    + sizeof(Node *)            // _coalescence_node
    + 4*sizeof(set<unsigned>)   // _left_species, _right_species, _anc_species, _species
    + sizeof(map<set<unsigned>, unsigned>)  // _lineage_counts
);
unsigned long                       Epoch::_max_bytes_any_epoch = Epoch::_bytes_base_epoch;
#endif
        
unsigned                            Forest::_nstates            = 4;

unsigned                            Forest::_ntaxa              = 0;
vector<string>                      Forest::_taxon_names;

unsigned                            Forest::_nspecies           = 0;
vector<string>                      Forest::_species_names;
map<unsigned,unsigned>              Forest::_nexus_taxon_map;

unsigned                            Forest::_ngenes             = 0;
vector<string>                      Forest::_gene_names;

map<string, unsigned>               Forest::_taxon_to_species;

double                              Forest::_theta              = 0.05;
double                              Forest::_lambda             = 1.0;

double                              Forest::_theta_prior_mean   = 0.05;
double                              Forest::_lambda_prior_mean  = 1.0;

bool                                Forest::_update_theta       = true;
bool                                Forest::_update_lambda      = true;

double                              Forest::_small_enough       = 0.00001;

static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 required in order to use infinity()");
double                              Forest::_infinity = numeric_limits<double>::infinity();
double                              Forest::_negative_infinity = -numeric_limits<double>::infinity();

vector<double>                      Forest::_cumprobs;

PartialStore::leaf_partials_t       GeneForest::_leaf_partials;

const double                        Node::_smallest_edge_length = 1.0e-12;

string                              Proj::_program_name         = "proj";
unsigned                            Proj::_major_version        = 1;
unsigned                            Proj::_minor_version        = 0;

// version 1.0: samples from species tree distribution conditional on supplied gene trees using SMC

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
#if defined(USING_SIGNPOSTS)
    log_handle  = os_log_create("edu.uconn.eeb.phylogeny", OS_LOG_CATEGORY_POINTS_OF_INTEREST);
    signpost_id = os_signpost_id_generate(log_handle);
    assert(signpost_id != OS_SIGNPOST_ID_INVALID);
#endif

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
        output(format("Exception: %s\n") % x.what(), 2);
        output("Aborted.\n", 2);
        normal_termination = false;
    }
    catch(...) {
        output("Exception of unknown type!\n", 2);
        normal_termination = false;
    }
    
#if defined(USING_MPI)
	MPI_Finalize();
#endif

#if defined(LOG_MEMORY)
    if (normal_termination) {
        proj.memoryReport(memfile);
        ps.memoryReport(memfile);
        Epoch::memoryReport(memfile);
    }
    memfile.close();
#endif

    return 0;
}
